#include "IOHandler.h"
#include "Split.h"

#include <algorithm>     // need transform
#include <iostream>
#include <iomanip>

// read run parameters and store in map
int IOHandler::readParamFile(char *filename)  {

  ifstream parin;
  parin.open(filename,ios::in);
  if ( !parin.is_open() )  {
    cerr << "error: can't open parameter file " << filename << endl;
    exit(1);
  }

  string line, linecpy;
  vector<string> elements;
  while ( !parin.eof() )  {

    getline(parin,line);
    linecpy = line;
    Split(line,&elements,":");

    if ( elements.size() == 0 )  {
      break;
    }
    else if ( elements.size() == 2 )  {
      string pname = elements[0];
      string pval = elements[1];

      map<string,string>::iterator p;
      p = parameters.find(pname);
      if ( p != parameters.end() )  {
        (*p).second = pval;
        cerr << "warning: redefinition of parameter " << pname << endl;
      }
      else  {
        parameters[pname] = pval;
      }
    }
    else  {
      cerr << "warning: ignoring unparsable parameter line " << endl 
           << linecpy << endl;
    }

  }
  parin.close();
  return 0;

}


int IOHandler::readIndFile()  {

  string fname = parameters["ind"];
  ifstream parin;
  parin.open(fname.c_str(),ios::in);
  if ( !parin.is_open() )  {
      cerr << "error: can't open indiv file " << fname << endl;
      exit(1);
  }

  vector<string> elements;
  nindiv = 0;
  while ( !parin.eof() )  {
    string line, linecpy;

    getline(parin,line);
    linecpy = line;               // tokenizer will mangle this but error handler may need it
    Split(line,&elements);

    if ( elements.size() == 0 )  {
      break;
    }
    else if ( elements.size() >= 2 )  {
      indivnames.push_back(elements[0]);
      indivpops.push_back(elements[2]);
      nindiv++;
    }
    else  {
      cerr << "warning: ignoring unparsable individual file line " << endl 
           << linecpy << endl;
    }
  }

  parin.close();
  return 0;

}

int IOHandler::readSnpWtFile()  {
    ifstream parin;
    string fname = parameters["snpwt"];
    parin.open(fname.c_str(),ios::in);
    if ( !parin.is_open() )  {
        cerr << "error: can't open snpwt file " << fname << endl;
        exit(1);
    }
    
    nwtsnps = 0;
    int linenum = 0;
    vector<string> elements;
    while ( !parin.eof() )  {

        string line, linecpy;
        
        getline(parin,line);
        linecpy = line;               // tokenizer will mangle this but error handler may need it
        Split(line,&elements);

        if ( elements.size() == 0 )  {
          break;
        }
    
        // first line is shrinkage one element per PC
        if ( linenum == 0 )  {
          npcs = elements.size();
          shrinkage = new double[npcs];
          for(int i=0;i<npcs;i++)  {
            shrinkage[i] = atof(elements[i].c_str());
          }
          linenum++;
          continue;
        }
        
        // second line is names of populations
        if ( linenum == 1 )  {
          npop = elements.size();
          linenum++;
          continue;
        }

        // third and fourth lines are not used
        if ( linenum == 2 || linenum == 3 )  { 
            linenum++;
            continue;
        }

        // fifth line is coefficients
        if ( linenum == 4 )  {   
          for(int i=0;i<(npcs+1)*(npcs+1);i++)  {
            coeff.push_back(atof(elements[i].c_str()));
          }
          linenum++;
          continue;
        }

        
        // SNP data lines
        string rsid = elements[0];
        wtref[rsid] = elements[1];
        wtvar[rsid] = elements[2];
        pmaf[rsid] = atof(elements[3].c_str());

        double *w = new double[npop];  // TODO : dealloc these
        for(int i=0;i<npcs;i++) {
            w[i] = atof(elements[4+i].c_str());
        }
        weights[rsid] = w;
        nwtsnps++;
    }
    
    /* compute flip */
    map<string,string> strandtrans;
    strandtrans["A"] = "T";
    strandtrans["T"] = "A";
    strandtrans["C"] = "G";
    strandtrans["G"] = "C";

    // TODO : error checking - is this SNP in {wtref, inref, wtvar, invar} ?
    // What should we do if it's not? Remove it?

    vector<string>::size_type i;
    for(i=0;i<rsids.size();i++)  {
        int ifl = -1;
        if ( weights.find(rsids[i]) != weights.end() )  {
            string r0 = wtref[rsids[i]];
            string r1 = inref[rsids[i]];
            string v0 = wtvar[rsids[i]];
            string v1 = invar[rsids[i]];
            
            if ( ( r0.compare(r1) == 0 && v0.compare(v1) == 0 ) ||
                 ( r0.compare(strandtrans[r1]) == 0 && v0.compare(strandtrans[v1]) == 0 ) )  {
                ifl = 0;
            }
            else if ( ( r0.compare(v1) == 0 && v0.compare(r1) == 0 ) ||
                      ( r0.compare(strandtrans[v1]) == 0 && v0.compare(strandtrans[r1]) == 0 ))  {
                ifl = 1;
            }
            else  {
                cout << "warning: unidentifiable flip (monomorphic?) - will not be used " 
                     << rsids[i] << " (" << r0 << " " << v0 << ") (" << r1 << " " << v1 << ")\n";
            }
        }
        flip[rsids[i]] = ifl;
    }
    
    
    return 0;
}

int IOHandler::readSnpFile()  {

    ifstream parin;
    string fname = parameters["snp"];
    parin.open(fname.c_str(),ios::in);
    if ( !parin.is_open() )  {
        cerr << "error: can't open indivfile " << fname << endl;
        exit(1);
    }

    vector<string> elements;
    ninsnps = 0;
    while ( !parin.eof() )  {

        string line, linecpy, rsid;

        getline(parin,line);
        linecpy = line;               // tokenizer will mangle this but error handler may need it
        Split(line,&elements);

        if ( elements.size() == 0 )  {
          break;
        }

        rsid = elements[0];
        rsids.push_back(rsid);
        inref[rsid] = elements[4];
        invar[rsid] = elements[5];
        ninsnps++;
    }

    parin.close();
    return 0;

}

int IOHandler::writeOutFile()  {

    ofstream outstream;
    string fname = parameters["predpcoutput"];
    outstream.open(fname.c_str());

    int npop = npcs+1;
    double *anc = new double[npop];

    if ( !outstream.is_open() )  {
        cerr << "error: can't open output file " << fname << endl;
        return 1;
    }

    for(int i=0;i<nindiv;i++)  {

        for(int j=0;j<npop;j++)
          anc[j] = 0.0;

        for(int j=0;j<npop;j++)  {
          for(int k=0;k<npcs;k++)  {
            anc[j] += coeff[j*npop + k]*Wout[i*npcs+k];
          }
          anc[j] += coeff[j*npop + npcs];
          if ( anc[j] < 0.0 ) anc[j] = 0.0;
          if ( anc[j] > 100.0 ) anc[j] = 100.0;
        }

        double sum = 0.0;
        for(int j=0;j<npop;j++)  {
          sum += anc[j];
        }
        if ( sum > 1.0 )  {
          for(int j=0;j<npop;j++)  {
            anc[j] /= sum;
          }
         }

        outstream << indivnames[i] << " " << indivpops[i] << " ";
        outstream << Nout[i] << " ";
        for(int k=0;k<npcs;k++)  {
           outstream << setprecision(10) << Wout[i*npcs+k] << " ";
        }
        for(int k=0;k<npop;k++)  {
           outstream << setprecision(10) << anc[k] << " ";
        }
        outstream << endl;
    }
    outstream.close();
    delete anc;
    return 0;
}


// Debugging stuff
int IOHandler::showParams()  {
    map<string,string>::iterator iter;
    for(iter=parameters.begin();iter != parameters.end();iter++)  {
        cout << iter->first << " " << iter->second << endl;
    }
    
    
    return 0;
}

int IOHandler::showIndivs()  {
    vector<string>::size_type i;
    for (i=0;i<=indivnames.size();i++)  {
        cout << indivnames[i] << " " << indivpops[i] << endl;
    }
    return 0;
}

int IOHandler::showSnps()  {
    vector<string>::size_type i;
    for(i=0;i<20;i++)  {
        cout << rsids[i] << " " << inref[rsids[i]] << " " << invar[rsids[i]] << endl;
    }
    return 0;
}

int IOHandler::showWeights()  {
    vector<string>::size_type i;
    for(i=0;i<40;i++)  {
        if ( weights.find(rsids[i]) != weights.end() )  {
            cout << rsids[i] << "  " << wtref[rsids[i]] << " " << wtvar[rsids[i]] << " "
                 << pmaf[rsids[i]] << "  ";
            for(int j=0;j<npcs;j++)  {
                cout << weights[rsids[i]][j] << " ";
            }
            cout << endl;
        }
        
    }
    return 0;
}









