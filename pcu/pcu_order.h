/****************************************************************************** 

  Copyright 2014 Scientific Computation Research Center, 
      Rensselaer Polytechnic Institute. All rights reserved.
  
  This work is open source software, licensed under the terms of the
  BSD license as described in the LICENSE file in the top-level directory.

*******************************************************************************/
#ifndef PCU_ORDER_H
#define PCU_ORDER_H

#include <stdbool.h>
#include "pcu_msg.h"

typedef struct pcu_order_struct* pcu_order;

pcu_order pcu_order_new(void);
void pcu_order_free(pcu_order o);
bool pcu_order_receive(pcu_order o, pcu_msg* m);
void* pcu_order_unpack(pcu_order o, size_t size);
bool pcu_order_unpacked(pcu_order o);
int pcu_order_received_from(pcu_order o);
size_t pcu_order_received_size(pcu_order o);

#endif
