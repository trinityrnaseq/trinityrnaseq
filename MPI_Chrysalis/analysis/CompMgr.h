#ifndef COMPMGR_H
#define COMPMGR_H

#include <string.h>

extern "C" {
#include <sys/types.h>
#include <sys/stat.h>
#include <unistd.h>
}


class ComponentFileMgr
{
 public:
  ComponentFileMgr(const string & base = ".") {
    m_base = base + "/";
    m_per = 2000;
    m_dir = m_base + "RawComps.";
  }

  string GetFileName(int index, const string & suffix, bool checkForDirectory = true) {
    int i = index / m_per;
    char tmp[1024];
    sprintf(tmp, "%s%d", m_dir.c_str(), i);
    if(checkForDirectory)
    {
       mkdir(tmp,0777);
    }
    string out = tmp;
    sprintf(tmp, "/comp%d", index);
    out += tmp;
    out += suffix;
    return out;
  }


 private:
  string m_base;
  string m_dir;
  int m_per;

};


#endif

