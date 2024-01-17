//
//  InferAnc.h
//  SNPweights
//
//  Created by user-zero on 10/26/13.
//  Copyright (c) 2013 Self. All rights reserved.
//

#ifndef _INFERANC_
#define _INFERANC_

#include <iostream>
#include <string>
#include <cstdlib>
#include <stdio.h>    // fgets in here

#include "IOHandler.h"

using namespace std;

class InferAnc  {
    
    public:
    
    InferAnc(IOHandler *p)  {
        ioh_p = p;
    }
    
    void Run(void);
    
    private:
    
    IOHandler *ioh_p;
    double *WX;
    int *ncount;
    
};




#endif
