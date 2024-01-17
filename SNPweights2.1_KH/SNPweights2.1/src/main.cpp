//
//  main.cpp
//  SNPweights
//
//  Created by user-zero on 10/26/13.
//  Copyright (c) 2013 Self. All rights reserved.
//

#include <iostream>
#include <string>
#include <cstring>
#include "IOHandler.h"
#include "InferAnc.h"

using namespace std;

int main(int argc, const char * argv[])  {
    
    if ( argc != 3 )  {
      cerr << "usage: inferanc -p <param-file>" << endl;
    }
    IOHandler *ioh = new IOHandler;
    char pname[256];
    strcpy(pname,argv[2]);
    
    ioh->readParamFile(pname);
    // ioh->showParams();
    ioh->readIndFile();
    // ioh->showIndivs();
    ioh->readSnpFile();
    // ioh->showSnps();
    ioh->readSnpWtFile();
    // ioh->showWeights();
    
    InferAnc *ia = new InferAnc(ioh);
    ia->Run();
}

