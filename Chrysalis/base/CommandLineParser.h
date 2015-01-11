/* -*- c-basic-offset: 4; indent-tabs-mode: nil -*- */

#include <string>
#include <typeinfo>
#include <sstream>
#include <iostream>

#include <map>
#include <set>
#include <stdlib.h>

#include "base/StringUtil.h"
#include "base/FileParser.h"

using namespace std;

template<typename argType>
class commandArg
{

 public:

    commandArg() {}

    commandArg(string name, string descrip) {
        mName = name;
        mDesc = descrip;
        mHasDefault = false;
        mHasValue = false;
    }

    commandArg(string name,string descrip,argType T) {
        mName = name;
        mDesc = descrip;
        mValue = T;
        mHasDefault = true;
        mHasValue = true;
    }

    void SetValue(argType &T) { mValue = T; mHasValue=true; }

    void SetDefault(argType &T) { mValue = T; mHasDefault = true; }
    void SetDescription(string d) { mDesc = d; }

    bool HasDefault() { return mHasDefault; }
    bool HasValue() { return mHasValue; }

    argType GetValue() { return mValue; }
    string GetName() const { return mName; }
    string GetDescription() { return mDesc; }

    friend bool operator< (const commandArg<argType> &lhs,
                           const commandArg<argType> &rhs)
    {
        if ( lhs.GetName() < rhs.GetName() )
            return true;

        return false;
    }

    string GetType()
    {
        string type = typeid(*this).name();
        if ( type.find("Ib") != string::npos )
            return "bool";
        else if ( type.find("Ii") != string::npos )
            return "int";
        else if ( type.find("Id") != string::npos )
            return "double";
        else if ( type.find("ISs") != string::npos )
            return "string";
        else if (type.find("IlE") != string::npos )
            return "long";
        else {
            std::cerr << "*** unknown type for [" << type << "]" << endl;
            return "unknown type";
        }
    }


 private:

    string mName, mDesc;
    argType mValue;
    bool mHasDefault, mHasValue;

};


class commandLineParser
{

 public:

    commandLineParser(int argc, char** argv, const string & description = "") {
        mArgc = argc;
        mArgv = argv;
        mDesc = description;
        mNumDefaults = 0;
        mNumArgs = 0;
    }

    void SetDescription(const string & s) {
        mDesc = s;
    }

 private:
    int mArgc;
    char** mArgv;
    stringstream mHelp;
    string mDesc;

    int mNumDefaults,mNumArgs;
    multimap<string,string> mNameVal;
    set<string> mArgNames, mCArgNames;
 public:


    template<typename T>
        void registerArg(commandArg<T> &arg)
        {
            if ( mHelp.str().empty() )
                mHelp << "\n";
            mHelp << string(arg.GetName()
                            + "<"
                            + arg.GetType()
                            + "> : "
                            + arg.GetDescription() );
            if ( arg.HasDefault() )
                {
                    mHelp << " (def=";
                    ostringstream o;
                    o << arg.GetValue() ;
                    mHelp << o.str() << ")";
                    ++mNumDefaults;
                }
            mHelp << "\n" << ends;

            mArgNames.insert(arg.GetName());
            ++mNumArgs;

        }

    template<typename T>
        void registerCompoundArg(commandArg<T> &arg)
        {
            if ( mHelp.str().empty() )
                mHelp << "\n";
            mHelp << string(arg.GetName()
                            + "<"
                            + arg.GetType()
                            + "> : "
                            + arg.GetDescription() );
            if ( arg.HasDefault() )
                {
                    mHelp << " (def=";
                    ostringstream o;
                    o << arg.GetValue() ;
                    mHelp << o.str() << ")";
                    ++mNumDefaults;
                }
            mHelp << "\n" << ends;

            mCArgNames.insert(arg.GetName());
            ++mNumArgs;
        }


    bool parse()
    {
        if ( mArgc == 1 && mNumArgs != mNumDefaults)
            {
                showHelp();
                exit(-1);
            }

        if ( requestHelp() )
            {
                showHelp();
                exit(-1);
            }

        if ( mArgNames.empty() )
            {
                showHelp();
                exit(-1);
            }

        int i(1);
        while ( i<mArgc )
            {
                //   for ( int i=1; i<mArgc; ++i )
                //   {


                string n(mArgv[i]);//, v(mArgv[i+1]);
                string v("");
                if ( i<mArgc-1)
                    v=string(mArgv[i+1]);

                bool is_compound(false);
                set<string>::iterator sIter = mArgNames.find(n);
                if ( sIter == mArgNames.end() )
                    {
                        sIter = mCArgNames.find(n);
                        if (sIter==mCArgNames.end())
                            {
                                if (n == "-print-command-line") {
                                    cout << "----------------------- Welcome to Trinity -------------------------------" << endl;
                                    cout << "This module was invoked via:" << endl;
                                    for (int k=0; k<mArgc; k++) {
                                        cout << mArgv[k] << " ";
                                    }
                                    cout << endl;
                                    cout << "----------------------- Welcome to Trinity -------------------------------" << endl << endl;
                                    if ( mArgc == 2 ) {
                                        showHelp();
                                        exit(-1);
                                    }
                                    i++;
                                    continue;
                                } else {

                                    cout << "\nInvalid command-line arg: " << n << endl;
                                    showHelp();
                                    exit(-1);
                                }
                            }
                        else
                            {
                                is_compound=true;
                            }
                    }

                string targ("-");
                StringParser sp; sp.SetLine(v);
                int num_v=sp.GetItemCount();

                if (!is_compound)
                    {
                        bool is_float(false);
                        if (num_v>1) is_float = sp.IsFloat(1);

                        if (ContainsAt(n,targ,0) && is_float)
                            {
                                mNameVal.insert(make_pair(n,v));
                                i+=2;
                            }
                        else if ( ContainsAt(n,targ,0) && !ContainsAt(v,targ,0) )
                            {
                                mNameVal.insert(make_pair(n,v));
                                i+=2;
                            }
                        else if ( ContainsAt(n,targ,0) && ContainsAt(v,targ,0) )
                            {
                                v="";
                                mNameVal.insert(make_pair(n,v));
                                i+=1;
                            }
                        else
                            {
                                i+=1;
                            }
                    }
                else
                    {
                        for (int k=0; k<num_v; ++k)
                            mNameVal.insert(make_pair(n,sp.AsString(k)));
                        i+=2;
                    }
            }
        return true;
    }

    template<typename T>
        string GetStringValueFor(commandArg<T> &keyArg)
        {
            string key(keyArg.GetName());
            map<string,string>::iterator mIter = mNameVal.find(key);
            if ( mIter == mNameVal.end() )
                {
                    ostringstream o;
                    o << keyArg.GetValue();

                    if ( keyArg.GetType() == "string" )
                        {
                            if ( !keyArg.HasDefault() )
                                {
                                    cout << "need to specify " << key << endl;
                                    exit(-1);
                                }
                            return (o.str());
                        }
                    return (string(""));
                }

            return (mIter->second);

        }

    template<typename T>
        vector<string> GetCompoundStringValuesFor(commandArg<T> &keyArg)
        {
            vector<string> dummy;

            string key(keyArg.GetName());
            pair<multimap<string,string>::iterator,multimap<string,string>::iterator>  mIter = mNameVal.equal_range(key);
            if ( mIter.first == mIter.second)
                {
                    ostringstream o;
                    o << keyArg.GetValue();

                    if ( keyArg.GetType() == "string" )
                        {
                            if ( !keyArg.HasDefault() )
                                {
                                    cout << "need to specify " << key << endl;
                                    exit(-1);
                                }

                            dummy.push_back(o.str());
                            return (dummy);
                        }
                    return (dummy);
                }

            for (; mIter.first != mIter.second; ++mIter.first)
                dummy.push_back(mIter.first->second);

            return dummy;

        }


    double GetDoubleValueFor(commandArg<double> &keyArg)
    {
        string key(keyArg.GetName());
        string val = GetStringValueFor(keyArg);

        if ( val.empty() )
            {
                if (!keyArg.HasDefault() )
                    {
                        cout << "need to specify " << key << endl;
                        exit(-1);
                    }
                else
                    return ( keyArg.GetValue() );
            }
        else
            {
                //     if ( !val.IsDouble() )
                //     {
                //       cout << "invalid: " << key <<" "<< val << endl;
                //       exit(-1);
                //     }
                //     else
                return (atof(val.c_str()));

            }
    }

    int GetIntValueFor(commandArg<int> &keyArg)
    {
        string key(keyArg.GetName());
        string val = GetStringValueFor(keyArg);
        if ( val.empty() )
            {
                if (!keyArg.HasDefault() )
                    {
                        cout << "need to specify " << key << endl;
                        exit(-1);
                    }
                else
                    return ( keyArg.GetValue() );
            }
        else
            {
                //     if ( !val.IsInt() )
                //     {
                //       cout << "invalid: " << key <<" "<< val << endl;
                //       exit(-1);
                //     }
                //     else
                return (atoi(val.c_str()));

            }
    }

    int GetLongValueFor(commandArg<long> &keyArg)
    {
        string key(keyArg.GetName());
        string val = GetStringValueFor(keyArg);
        if ( val.empty() )
            {
                if (!keyArg.HasDefault() )
                    {
                        cout << "need to specify " << key << endl;
                        exit(-1);
                    }
                else
                    return ( keyArg.GetValue() );
            }
        else
            {
                //     if ( !val.IsInt() )
                //     {
                //       cout << "invalid: " << key <<" "<< val << endl;
                //       exit(-1);
                //     }
                //     else
                return (atol(val.c_str()));

            }
    }







    bool GetBoolValueFor(commandArg<bool> &keyArg)
    {
        string key(keyArg.GetName());
        map<string,string>::iterator mIter = mNameVal.find(key);
        if ( mIter == mNameVal.end() )
            {
                if ( keyArg.HasDefault() )
                    {
                        return keyArg.GetValue();
                    }
                else
                    {
                        cout << "need to specify " << key << endl;
                        exit(-1);
                    }
            }
        else
            return true;
    }


    bool requestHelp()
    {
        for ( int i=1; i<mArgc; ++i )
            {
                string this_arg(mArgv[i]);
                if ( this_arg == "-h" )
                    return true;
            }

        return false;
    }

    void showHelp()
    {
        cout << endl << mArgv[0] << ": ";
        if (mDesc != "")
            cout << mDesc << endl << endl;
        else
            cout << "a module in the code base 'Trinity'." << endl << endl;
        cout << "\nAvailable arguments:" << endl;
        cout << mHelp.str() << endl;
    }


};

