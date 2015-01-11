#include "string_util.hpp"

namespace string_util {
    
    void tokenize(const string& str, vector<string>& tokens, const string& delimiters) {
        
        /* ************************************************************************************/
        /* from: http://oopweb.com/CPP/Documents/CPPHOWTO/Volume/C++Programming-HOWTO-7.html **/
        /**************************************************************************************/
        
        // Skip delimiters at beginning.
        string::size_type lastPos = str.find_first_not_of(delimiters, 0);
        // Find first "non-delimiter".
        string::size_type pos     = str.find_first_of(delimiters, lastPos);
        
        while (string::npos != pos || string::npos != lastPos) {
            // Found a token, add it to the vector.
            tokens.push_back(str.substr(lastPos, pos - lastPos));
            // Skip delimiters.  Note the "not_of"
            lastPos = str.find_first_not_of(delimiters, pos);
            // Find next "non-delimiter"
            pos = str.find_first_of(delimiters, lastPos);
        }
        
        return;
    }
    

    /* rest is my stuff */

    string join(const vector<string>& tokens, const string& delimiter) {
        
        string concat = "";
        
        for (unsigned int i = 0; i < tokens.size(); i++) {
            concat = concat + tokens[i];
            
            if (i != tokens.size()-1) {
                concat = concat + delimiter;
            }
        }

        return(concat);
    }
    




} // end of string_util namespace

