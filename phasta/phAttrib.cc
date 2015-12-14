# include "phAttrib.h"
#include "gmi_sim.h"
#include <SimAttribute.h>
#include <SimUtil.h>
#include <cstdlib>
#include <iostream>

/* Simmetrix, for the love of all that is good,
   please put this in your header files.
   We can't work with files from Simmodeler without it. */
pAManager SModel_attManager(pModel model);

typedef ph::BC* (*BCFactory)(pAttribute a, pGEntity ge);
typedef std::map<std::string, BCFactory> BCFactories;

struct SimBC : public ph::BC
{
  SimBC(pGEntity ge)
  {
    dim = GEN_type(ge);
    tag = GEN_tag(ge);
  }
};

struct Tensor0BC : public SimBC
{
  Tensor0BC(pAttribute a, pGEntity ge):SimBC(ge)
  {
    if (Attribute_repType(a) != Att_tensor0) {
      fprintf(stderr, "tensor 0 attribute does not match type\n");
      abort();
    }
    attribute = (pAttributeTensor0)a;
  }
  virtual double* eval(apf::Vector3 const& x)
  {
    buf = AttributeTensor0_evalDS(attribute, &x[0]);
    return &buf;
  }
  pAttributeTensor0 attribute;
  double buf;
};

struct Tensor1BC : public SimBC {
  Tensor1BC(pAttribute a, pGEntity ge):SimBC(ge)
  {
    if (Attribute_repType(a) != Att_tensor1) {
      fprintf(stderr, "tensor 1 attribute does not match type\n");
      abort();
    }
    attribute = (pAttributeTensor1)a;
  }
  virtual double* eval(apf::Vector3 const& x)
  {
    for (int i = 0; i < 3; ++i)
      buf[i] = AttributeTensor1_evalDS(attribute, i, &x[0]);
    return buf;
  }
  pAttributeTensor1 attribute;
  double buf[3];
};

struct CompBC : public SimBC {
  CompBC(pAttribute a, pGEntity ge):SimBC(ge)
  {
    if (Attribute_repType(a) != Att_void) {
      fprintf(stderr, "comp1/3 attribute does not match type\n");
      abort();
    }
    pPList children = Attribute_children(a);
    magnitude = 0;
    direction = 0;
    for (int i = 0; i < PList_size(children); ++i) {
      pAttribute child = (pAttribute) PList_item(children, i);
      if (Attribute_repType(child) == Att_double)
        magnitude = (pAttributeDouble) child;
      else if (Attribute_repType(child) == Att_tensor1)
        direction = (pAttributeTensor1) child;
else
  fprintf(stderr,"ignored some comp1/3 attributes...\n");
    }
    PList_delete(children);
    if (!magnitude) {
      fprintf(stderr, "comp1/3 attribute does not have magnitude\n");
      abort();
    }
    if (!direction) {
      fprintf(stderr, "comp1/3 attribute does not have direction\n");
      abort();
    }
  }
  virtual double* eval(apf::Vector3 const& x)
  {
    buf[0] = AttributeDouble_evalDS(magnitude, &x[0]);
    for (int i = 0; i < 3; ++i)
      buf[i + 1] = AttributeTensor1_evalDS(direction, i, &x[0]);
    return buf;
  }
  pAttributeDouble magnitude;
  pAttributeTensor1 direction;
  double buf[4];
};

struct IntBC : public SimBC
{
  IntBC(pAttribute a, pGEntity ge):SimBC(ge)
  {
    if (Attribute_repType(a) != Att_int) {
      fprintf(stderr, "int attribute does not match type\n");
      abort();
    }
    attribute = (pAttributeInt)a;
  }
  virtual double* eval(apf::Vector3 const&)
  {
    /* forgive me for I am storing an int in a double.
       let there be more than 32 mantissa bits such that
       this conversion is lossless. */
    buf = AttributeInt_value(attribute);
    return &buf;
  }
  pAttributeInt attribute;
  double buf;
};

static ph::BC* tensor0Factory(pAttribute a, pGEntity ge)
{
  return new Tensor0BC(a, ge);
}

static ph::BC* tensor1Factory(pAttribute a, pGEntity ge)
{
  return new Tensor1BC(a, ge);
}

static ph::BC* compFactory(pAttribute a, pGEntity ge)
{
  return new CompBC(a, ge);
}

static ph::BC* intFactory(pAttribute a, pGEntity ge)
{
  return new IntBC(a, ge);
}

/* this should follow the KnownBC tables in phBC.cc */
static void formFactories(BCFactories& fs)
{
  fs["density"]              = tensor0Factory;
  fs["temperature"]          = tensor0Factory;
  fs["pressure"]             = tensor0Factory;
  fs["comp1"]                = compFactory;
  fs["comp3"]                = compFactory;
  fs["scalar_1"]             = tensor0Factory;
  fs["scalar_2"]             = tensor0Factory;
  fs["scalar_3"]             = tensor0Factory;
  fs["scalar_4"]             = tensor0Factory;
  fs["mass flux"]            = tensor0Factory;
  fs["natural pressure"]     = tensor0Factory;
  fs["traction vector"]      = tensor1Factory;
  fs["heat flux"]            = tensor0Factory;
  fs["turbulence wall"]      = tensor0Factory;
  fs["scalar_1 flux"]        = tensor0Factory;
  fs["scalar_2 flux"]        = tensor0Factory;
  fs["scalar_3 flux"]        = tensor0Factory;
  fs["scalar_4 flux"]        = tensor0Factory;
  fs["surf ID"]              = tensor0Factory;
  fs["initial pressure"]     = tensor0Factory;
  fs["initial velocity"]     = tensor1Factory;
  fs["initial temperature"]  = tensor0Factory;
  fs["initial scalar_1"]     = tensor0Factory;
  fs["initial scalar_2"]     = tensor0Factory;
  fs["initial scalar_3"]     = tensor0Factory;
  fs["initial scalar_4"]     = tensor0Factory;
  fs["periodic slave"]       = intFactory;
  fs["DG interface"]         = intFactory;
  fs["material type"]        = intFactory;
}

static void addAttribute(BCFactories& fs, pAttribute a, pGEntity ge,
    ph::BCs& bcs)
{
  char* c_infoType = Attribute_infoType(a);
  std::string infoType(c_infoType);
  if (!fs.count(infoType)) {
    fprintf(stderr,"unknown attribute type \"%s\", ignoring !\n", c_infoType);
    fprintf(stderr,"it had repType %d\n",
        Attribute_repType(a));
    return;
  }
  if (!bcs.fields.count(infoType))
    bcs.fields[infoType] = ph::FieldBCs();
  ph::FieldBCs& fbcs = bcs.fields[infoType];
  Sim_deleteString(c_infoType);
  BCFactory f = fs[infoType];
  fbcs.bcs.insert( f(a, ge) );
}

static void addAttributes(BCFactories& fs, pPList as, pGEntity ge,
    ph::BCs& bcs)
{
  for (int i = 0; i < PList_size(as); ++i) {
    pAttribute a = (pAttribute) PList_item(as, i);
    addAttribute(fs, a, ge, bcs);
  }
}

void ph::getSimmetrixAttributes(gmi_model* model, ph::BCs& bcs)
{
  pGModel smdl = gmi_export_sim(model);
  pAManager mngr = SModel_attManager(smdl);
  if (!mngr) {
    fprintf(stderr,"Simmetrix model did not come with an Attribute Manager\n");
    abort();
  }
  pACase pd = AMAN_findCaseByType(mngr, "problem definition");
  if (!pd) {
    fprintf(stderr,"no Attribute Case \"problem definition\"\n");
    abort();
  }
  AttCase_setModel(pd, smdl);
  AttCase_associate(pd, NULL);
  BCFactories fs;
  formFactories(fs);
  for (int dim = 0; dim < 4; ++dim) {
    gmi_iter* it = gmi_begin(model, dim);
    gmi_ent* e;
    while ((e = gmi_next(model, it))) {
      pGEntity ge = (pGEntity) e;
      pPList attribs = GEN_attributes(ge, "");
      addAttributes(fs, attribs, ge, bcs);
    }
    gmi_end(model, it);
  }
  /* Simmodeler seems to put ICs on the model itself, not the regions.
     spread all ICs off the model onto the regions.
     also, this code does assume the domain is 3D. */
  pPList attribs = GM_attributes(smdl, "");
  GRIter rit = GM_regionIter(smdl);
  pGRegion gr;
  while ((gr = GRIter_next(rit)))
    addAttributes(fs, attribs, gr, bcs);
  GRIter_delete(rit);
}

