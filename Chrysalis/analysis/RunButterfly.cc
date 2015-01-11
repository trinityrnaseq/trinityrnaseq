#include <string>
#include <unistd.h>
#include "base/CommandLineParser.h"
#include "base/FileParser.h"


class JobControl
{
public:
  JobControl(int n, bool bLSF) {
    m_n = n;
    m_run = 0;
    m_bLSF = bLSF;
  }
  
  void Add(const string & cmmd, const string & out, const string & log, const string &err) {
    m_cmmd.push_back(cmmd);
    m_out.push_back(out);
    m_log.push_back(log);
    m_err.push_back(err);
    m_sub.push_back(0);
  }


  int Remaining() {
    return m_cmmd.isize();      
  }

  void Submit() {
    int i;

    // Check if they are done...
    for (i=0; i<m_sub.isize(); i++) {
      if (m_sub[i] == 0)
	continue;
      FILE * p = fopen(m_err[i].c_str(), "r");
      if (p == NULL)
	continue;
      
      fclose(p);

      FlatFileParser parser;
  
      parser.Open(m_err[i]);

      bool bDone = false;
      while (parser.ParseLine()) {
	if (parser.Line() == "Done") {
	  bDone = true;
	}
      }

      if (bDone) {
	Remove(i);
	cout << "JC: Done, removing: " << i << endl;
	m_run--;
	i--;
      }

    }
    //cout << "JC: Running: " << m_run << endl;
    //cout << "JC: Possibly submitting " << m_n - m_run << " new jobs." << endl;
  
    if (m_run == m_n)
      return;

    for (i=0; i<m_sub.isize(); i++) {
      if (m_sub[i] == 0) {
	m_sub[i] = 1;
	m_run++;
	
	if (!m_bLSF) {
	  string theCmmd = m_cmmd[i];
	  theCmmd += " >& ";
	  theCmmd += m_log[i];
	  theCmmd += " &";
	  system(theCmmd.c_str());
	  cout << "JC: Start job: " << theCmmd << endl;
	} else {
	  string theCmmd = "bsub -o ";
	  theCmmd += m_log[i];
	  theCmmd += " ";
	  theCmmd += m_cmmd[i];
	  cout << "JC: Start job on farm: " << theCmmd << endl;
	  system(theCmmd.c_str());	  
	}

	if (m_run >= m_n)
	  break;
      }
    }
    
  }

  void ClearFiles() {

    for (int i=0; i<m_err.isize(); i++) {
      //if (m_sub[i] == 0)
      //continue;
      FILE * p = fopen(m_err[i].c_str(), "r");
      if (p == NULL)
	continue;
      fclose(p);
      string rm = "rm " + m_err[i];
      cout << "Running: " << rm << endl;
      cout << rm << endl;
      system(rm.c_str());
    }

    
  }


  void Remove(int i) {
   
    int n = m_cmmd.isize()-1;
    if (n < 0)
      return;

    m_cmmd[i] = m_cmmd[n];
    m_out[i] = m_out[n];
    m_sub[i] = m_sub[n];
    m_log[i] = m_log[n];
    m_err[i] = m_err[n];
    m_cmmd.resize(n);
    m_sub.resize(n);
    m_out.resize(n);
    m_log.resize(n);
    m_err.resize(n);
  }

private:
  bool m_bLSF;
  int m_n;
  int m_run;
 
  svec<int> m_sub;
  svec<string> m_cmmd;
  svec<string> m_out;
  svec<string> m_log;
  svec<string> m_err;
};



int main( int argc, char** argv )
{


  commandArg<string> fileCmmd("-i","input file");
  commandArg<int> nCmmd("-n","number of CPUs", 1);
  commandArg<bool> lsfCmmd("-lsf","run on LSF or not",false);
  commandLineParser P(argc,argv);
  P.SetDescription("Parallel execution of butterfly on a local server or on the LSF farm.");
  P.registerArg(fileCmmd);
  P.registerArg(nCmmd);
  P.registerArg(lsfCmmd);
  
  P.parse();
  
  string fileName = P.GetStringValueFor(fileCmmd);
  int nCPUs = P.GetIntValueFor(nCmmd);
  bool bLSF = P.GetBoolValueFor(lsfCmmd);
  

  //comment. ???
  FlatFileParser parser;
  
  parser.Open(fileName);

  JobControl jc(nCPUs, bLSF);
  jc.ClearFiles();


  cout << "Reading command list..." << endl;




  //!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  //parser.ParseLine();
  while (parser.ParseLine()) {
    if (parser.GetItemCount() == 0)
      continue;
    string out = parser.AsString(parser.GetItemCount()-1);
    string log = out;
    string err = out;
    out += "_allProbPaths.fasta";
    log += ".log";
    err += ".err";
    jc.Add(parser.Line(), out, log, err);
  }  

  jc.ClearFiles();

  //return 0;


  do {
    //cout << "Submitting a batch." << endl;
    jc.Submit();
    //usleep(100000);
    sleep(1);
  } while (jc.Remaining() > 0);


  return 0;
}
