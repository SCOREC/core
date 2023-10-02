//============================================================
// This file is part of the Meshing and Geometry (MG) Project
// of the Computational Research and Engineering Acquisition
// Tools and Environments (CREATE) program.
//
// Use and/or redistribution is governed by agreement on file 
// with the CREATE Program Office. 
// Any unauthorized use and/or redistribution is prohibited.
//============================================================
/*! \file Main.cpp
 *  
 *  Test example of the use of Capstone Module back-end.
 *  This is a black box, no-GUI approach to using capstone tools
 *
 */

/** \example testMeshImplant.cpp
 Demonstrates mesh implant and remeshing for quality control using CapstoneModule.
 */

#include "CapstoneModule.h"
#include "CreateMG_Framework_Core.h"

#include "CreateMG_AttributionDatabaseInterface.h"
#include "CreateMG_GeometryDatabaseInterface.h"
#include "CreateMG_GeometryUtility.h"
#include "CreateMG_GeometrySmartIterator.h"
#include "CreateMG_DataDir.h"
#include "CreateMG_Error.h"
#include "CreateMG_Function.h"
#include "CreateMG_MeshDatabaseInterface.h"
#include "CreateMG_AttributesTool.h"
#include "CreateMG_Analysis.h"

using namespace CreateMG;
using namespace CreateMG::Geometry;
using namespace CreateMG::Mesh;
using namespace CreateMG::Attribution;
using namespace CreateMG::Utility;


MStatus demo_discrete_boolean(CapstoneModule &cs)
{
    GeometryDatabaseInterface *gdbi = cs.get_geometry();
    MeshDatabaseInterface *mdbi = cs.get_mesh();
    AttributionDatabaseInterface *adbi = cs.get_attribution();
    
    AppContext *ctx = cs.get_context();

    try
    {
        // 0. Create a mesh model and a geometry model
        M_GModel gmodel;
        M_MModel mmodel;
        //MG_API_CALL(gdbi, create_model(gmodel));
        //MG_API_CALL(mdbi, create_associated_model(mmodel, gmodel));
        //MG_API_CALL(mdbi, set_support_mesh_model(mmodel));
        //MG_API_CALL(gdbi, set_current_model(gmodel));
        //MG_API_CALL(mdbi, set_current_model(mmodel));
        
        // 0.1 Extract gmodel name to insert into the reader
        //std::string modelName;
        //MG_CALL(gdbi->get_model_name(gmodel, modelName));   
        
        // 1. Load two mesh files for terrain and component breps from stl or vtk files

        std::vector<std::string> files(2);
        files[0] = DataDir::get_data_file("chamfer_refined.stl");
        files[1] = DataDir::get_data_file("cube_2.stl");
        gmodel = cs.load_files(files);
        MG_API_CALL(gdbi, set_model_name(gmodel, "ImplantTest"));
        
        std::vector<M_MModel> mmodels;
        MG_API_CALL(mdbi, get_associated_mesh_models(gmodel, mmodels));
        
        if(mmodels.empty())
        {
            AppBBoard::error("No mesh models found");
            return ERR_IO;
        }

        mmodel = mmodels[0];
        MG_API_CALL(mdbi,set_current_model(mmodel));

        // 2. Check if we have two breps
        int num_breps;
        MG_API_CALL(gdbi, get_num_breps(num_breps));
        if(num_breps != 2)
        {
            AppBBoard::error("Does not have 2 breps");
            return ERR_IO;
        }
        
        M_GBRep brep[2];
        MG_API_CALL(gdbi, get_brep_by_index(0, brep[0]));
        MG_API_CALL(gdbi, get_brep_by_index(1, brep[1]));
        
        // 2. Translate box to intersect with the upper plane of the chamfer
        TransformationMatrix trans;
        trans.translation(2, 2, 1.5);
        MG_API_CALL(gdbi, apply_transform(brep[1], trans));
        
        // ALTERNATIVE WAY TO 3 BELOW THAT RUNS GENERIC BOOLEAN FUNCTION 
        // INSTEAD OF THE DISCRETE IMPLANT FUNCTION ; 
        // FIX: Issue found - gmodel and mmodel is retrieved from ui state which is not existing in this mode.
        /* 
         FunctionPtr boolean = get_function(context, "Boolean");
         if(!boolean.get())
         {
         AppBBoard::error("Could not find function boolean");
         return ERR_GENERIC;
         }
         // union, nonreg_union, difference, intersection, imprint
         set_input_model(boolean, gmodel);
         set_input(boolean,"Type", "nonreg_union");
         set_input(boolean,"Brep1", brep[0]);
         set_input(boolean,"Brep2", brep[1]);
         set_input(boolean,"KeepBreps", false);
         
         boolean->execute(context);
         
         M_GBRep obrep;
         get_output(boolean,"Output", obrep);        
         */
        
        
        // 3. Implant box (component-cube_2) into the terrain (ship-chamfer) brep 
        AppBBoard::highlight("Running Discrete Boolean from mesh implant function");
        FunctionPtr fct = cs.get_function("MeshImplant");
        if(!fct.get())
        {
            AppBBoard::error("No implant function.");
            return ERR_GENERIC;
        }
        
        // 3.1. set nonreg_union option
        std::string _type = "nonreg_union";
        int boolcode = -1;
        if(_type == "implant"     ) boolcode = 0;//IMPLANT (NONREGUNION-IMPRINT)
        if(_type == "nonreg_union") boolcode = 5;//NONREGUNION;
        if(_type == "union"       ) boolcode = 1;//UNION;
        if(_type == "difference"  ) boolcode = 2;//DIFFERENCE;
        if(_type == "intersection") boolcode = 3;//INTERSECTION;
        
        if(boolcode == -1)
        {
            AppBBoard::warning("Wrong boolean code");
            return ERR_INVALID_INPUT;
        }
        
        set_input(fct, "MeshModel"     ,mmodel);
        set_input(fct, "ComponentBrep" ,brep[1]);
        set_input(fct, "ShipBrep"      ,brep[0]);
        set_input(fct, "ImproveQuality",false);
        set_input(fct, "BooleanCode"   ,boolcode);
        
        cs.execute(fct);
        
        // 4. Create sizing scenario for remeshing the implanted brep
        
        M_GBRep obrep;
        get_output(fct,"NewBrep", obrep);

        //MG_API_CALL(gdbi, get_brep_by_index(0, obrep));
        
        if(obrep.is_invalid())
        {
            AppBBoard::error("There is no implanted brep found");
            return ERR_GENERIC;
        }
        
        // 4.1 Set a uniform size as a factor of brep's bounding box
        double size = 0.1;
        try {
            BBox bbox;
            GeometryUtility gt(gdbi);
            gt.get_bounding_box(std::vector<M_GBRep>(1, obrep), bbox);
            size = bbox.big_length() / 10.0;
        }
        catch(Exception &) {
        }
        
        Attribution::Analysis *analysis = get_analysis(ctx, "DefaultAnalysis");
        Attribution::AttributesTool tool(adbi, gmodel, analysis);
        tool.set_global_sizing(PD_GlobalSize, size);

        cs.set_meshers("MESHER C1");
        
        std::vector<M_GBRep> breps;
        breps.push_back(obrep);
        
        // 5. Execute a mesher with the sizing over the implanted brep
        FunctionPtr mesher = cs.get_function("GenerateFaceMesh");
        if(!mesher.get())
        {
            AppBBoard::error("No generate mesh function.");
            return ERR_GENERIC;
        }
        
        set_input_model(mesher, gmodel);
        set_input_analysis(mesher, analysis->get_name());
        set_input(mesher,"Breps" , breps);
        set_input(mesher,"InputMesh" , mmodel);
        
        cs.execute(mesher);
        
        std::vector<M_GTopo> failed;
        mesher->output().get_value("Failed", failed);
        if(!failed.empty())
        {
            AppBBoard::error("Failed to mesh %d faces.", failed.size());
            return ERR_GENERIC;
        }

        M_MModel ommodel;
        get_output(mesher, "MeshModel", ommodel);

        cs.save_mesh_file(SystemTool::get_real_app_path().append("/booleaned.vtk"),ommodel); 	

    }
    
    catch(CreateMG::Exception &err)
    {
        AppBBoard::error(err.output(ExceptionFormatter()).c_str());
        std::cout << "Implant Mesh failed." << std::endl;
        return err.get_status();
    }
    
    catch(...) {
        AppBBoard::error("Unknown error");
        return ERR_GENERIC;
    }

    /*
    std::string cre_1_0 =  SystemTool::get_real_app_path().append("/run_1_0.cre");
    if( cs.save_file(cre_1_0,gmodel1_0)!=STATUS_OK) {
        AppBBoard::warning("Could not save model %s",cre_1_0.c_str());
    }
    */

    return STATUS_OK;
}

/** Demonstrates mesh implant and remeshing for quality control using CapstoneModule.
 */
int main(int argc, char** argv)
{
    // extract the real location of the application executable (needed to find the plugins)
    SystemTool::set_real_app_path(argv[0]);
    
    // create a bulletin board for print-outs
    AppBBoard::add(new BBoardPosterConsole("console"));
    
    // set-up MG's exception mechanism
    // MG_CALL macro throws CreateMG::Exception if below is set to true
    // otherwise, if false, it does not throw exception, one can then list errors by Error::record_print(); 
    Error::exception(true);
    
    // create an instance of the Capstone Module activating SMLIB/CREATE/CREATE for the Geometry/Mesh/Attribution databases 
    const std::string gdbName("Geometry Database : Create");
    const std::string mdbName("Mesh Database : Create");
    const std::string adbName("Attribution Database : Create");
    CapstoneModule  cs("theModule", gdbName.c_str(), mdbName.c_str(), adbName.c_str());
    if (cs.get_geometry()==0 || cs.get_mesh()==0 || cs.get_attribution()==0) {
        AppBBoard::error("CapstoneModule could not be properly activated");
        if (cs.get_geometry()==0)
            AppBBoard::error("  Geometry : %s does not exist!",gdbName.c_str());
        if (cs.get_mesh()==0) 
            AppBBoard::error("  Mesh     : %s does not exist!",mdbName.c_str());
        if (cs.get_attribution()==0) 
            AppBBoard::error("  Attribution : %s does not exist!",adbName.c_str());
    }
    
    int result = 0;
    try {
        AppBBoard::normal("Discrete Booleans and Remeshing");
        AppBBoard::normal("-------------");
        if(demo_discrete_boolean(cs) != STATUS_OK)
            result = 1;
    }
    catch(Exception&) {
        AppBBoard::error("Error detected while executing Discrete Booleans and Remeshing");
        return 1;
    }
    return result;
}
