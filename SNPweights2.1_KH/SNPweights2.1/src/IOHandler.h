#ifndef _IOHANDLER_
#define _IOHANDLER_

#include <fstream>
#include <sstream>
#include <iomanip>
#include <cstdlib>

#include <vector>
#include <map>

using namespace std;

class IOHandler  {

    public:

    int readParamFile(char *filename);
    int readSnpFile();
    int readSnpWtFile();
    int readIndFile();
    int writeOutFile();
    
    int showParams();
    int showIndivs();
    int showSnps();
    int showWeights();
    
    
    // destructor - clean up allocated weight arrays
    
    // accessors
    int getNinSNPs()  {  return ninsnps;  }
    int getNwtSNPs()  {  return nwtsnps;  }

    int getNpop()     {  return npop;     }
    int getNpcs()     {  return npcs;     }
    int getNindiv()   {  return nindiv;   }
    
    double **getShrinkage()    {  return &shrinkage;  }

    map<string, double *> *getWeights()  {  return &weights;  }
    map<string, double> *getPmaf()       {  return &pmaf;     }
    map<string,int> *getFlip()           {  return &flip;     }
    vector<string> *getRsids()           {  return &rsids;    }
    
    string getGenoFileName()  {  return parameters["geno"];  }
    
    void setWout(double *p)  {  
        Wout = p;  
    }
    void setNout(int *p)     {  
        Nout = p;  
    }
    
    private:
    
    int ninsnps;
    int nwtsnps;
    int nindiv;
    int npcs;
    int npop;

    double* shrinkage;
    map<string,int> flip;

    vector<string> indivnames;
    vector<string> indivpops;
    vector<string> rsids;
    vector<double> coeff;

    map<string,string> inref;
    map<string,string> invar;
    map<string,string> wtref;
    map<string,string> wtvar;

    map<string,double *> weights;
    map<string,double> pmaf;

    map<string,string> parameters;
    
    double *Wout;
    int *Nout;


};

#endif
