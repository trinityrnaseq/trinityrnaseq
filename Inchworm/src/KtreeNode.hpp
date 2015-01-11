#ifndef __KtreeNode__

#define __KtreeNode__

#include <string>
#include <vector>
#include <map>

using namespace std;

class KtreeNode {

public:

    KtreeNode(char nodeChar, long id);
    
    long get_child(char nodeChar);

    vector<char> get_children();

    bool has_child(char nodeChar);
    
    void add_child(char nodeChar, long id);
        
    bool has_children();

    char get_char();
    
    long get_count();
    
    void set_count(long count);

    string toString();
    
    
private:

    char nodeChar;
        
    long count;
    
    long children[4]; // G, A, T, C => ID

    int char_to_position(char nodeChar); // G-> 0, A-> 1, ...

    char position_to_char(int position); // 0 -> G, 1-> A, ...

};

#endif
