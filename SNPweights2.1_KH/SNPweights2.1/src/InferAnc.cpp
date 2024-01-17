//
//  InferAnc.cpp
//  SNPweights
//
//  Created by user-zero on 10/26/13.
//  Copyright (c) 2013 Self. All rights reserved.
//

#include "InferAnc.h"
#include <cmath>

void InferAnc::Run()  {
    
    /* get pointers to data from IOHandler */
    int nindiv = ioh_p->getNindiv();
    int npc = ioh_p->getNpcs();
    int nsnps = ioh_p->getNinSNPs();
    
    double **shrinkage_p = ioh_p->getShrinkage();
    map<string, double *> *weights_p = ioh_p->getWeights();
    map<string,double> *pmaf_p = ioh_p->getPmaf();
    vector<string> *rsids_p = ioh_p->getRsids();
    map<string,int> *flip_p = ioh_p->getFlip();
    
    string genofilename = ioh_p->getGenoFileName();
    
    FILE *fid = fopen(genofilename.c_str(),"r");
    if ( fid == NULL )  {
        printf("error: can't open geno file %s\n", genofilename.c_str());
        exit(1);
    }
    
    char *buffer = new char[nindiv+2];    /* leave room for newline char and \0 */
    WX = new double[nindiv*npc];
    ncount = new int[nindiv];
    for(int j=0;j<nindiv*npc;j++)  {
        WX[j] = 0.0;
    }
    for(int j=0;j<nindiv;j++)  {
        ncount[j] = 0;
    }

    int idbg = 0;
    for(int j=0;j<nsnps;j++)  {
        
        string snp = (*rsids_p)[j];

        // if SNP is not in .snpwt file, read past the line 
        if ( (*weights_p).find(snp) == (*weights_p).end() )  {
            if ( fgets(buffer,nindiv+2,fid) != NULL )
                continue; 
            else
                break;     // EOF (should not be possible)
        }
        
        int iflip = (*flip_p)[snp];
        double *w = (*weights_p)[snp];
        double pmaf = (*pmaf_p)[snp];
        double factor = 1.0/sqrt(pmaf*(1.0-pmaf));

        if ( fgets(buffer,nindiv+2,fid) != NULL && iflip != -1 )  {

            for(int i=0;i<nindiv;i++)  {
                int ig = 9;
                if ( buffer[i] == '0' ) ig = 0;
                if ( buffer[i] == '1' ) ig = 1;
                if ( buffer[i] == '2' ) ig = 2;
                if ( ig != 9 && iflip == 1 )
                    ig = 2 - ig;
                if ( ig != 9 )  {
                    for(int k=0;k<npc;k++)  {
                        WX[i*npc+k] += w[k]*((double)ig-2.0*pmaf)*factor;
                    }
                    ncount[i] += 1;
                }
            }
        }
    }

    int nw = ioh_p->getNwtSNPs();
    for(int i=0;i<nindiv;i++)  {
        for(int k=0;k<npc;k++)  {
            WX[i*npc+k] *= ((double)nw/(double)ncount[i]);
            WX[i*npc+k] /= (*shrinkage_p)[k];
        }
    }

    ioh_p->setWout(WX);
    ioh_p->setNout(ncount);
    ioh_p->writeOutFile();
    
    return;

}






