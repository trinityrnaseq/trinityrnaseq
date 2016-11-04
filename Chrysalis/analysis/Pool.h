#ifndef __POOL_HPP__

#define __POOL_HPP__

// wrapper around an integer vector
class Pool
{
public:


    Pool() {
        m_id = -1;
    }

    Pool(int id) {
        m_id = id;
    }

    Pool(const Pool& p) {
        m_id = p.m_id;
        m_index = p.m_index;
    }

    void push_back(int i) { add(i); }
    void add(int i) {m_index.push_back(i);}

    void add (Pool& p) {
        for (int i = 0; i < p.size(); i++) {
            if (! contains(p[i])) {
                add(p[i]);
            }
        }
    }

    void clear () {
        svec<int> newvec;
        m_index = newvec;
    }

    void exclude (Pool& p) {

        for (int i = 0; i < p.size(); i++) {
            if (p.contains(i)) {
                exclude(i);
            }
        }
    }

    void exclude (int i) {
        //FIXME:  there must be a more efficient way to do this.
        svec<int> new_mvec;
        for (size_t j = 0; j < m_index.size(); j++) {
            int val = m_index[j];
            if (val != i) {
                new_mvec.push_back(val);
            }
        }
        m_index = new_mvec;
    }


    int size() const {return m_index.isize();}
    int isize() const { return(size()); }

    int get(int i) const {return m_index[i];}

    int get_id() {
        return(m_id);
    }

    int operator [] (const int& i) {
        return(get(i));
    }

    bool contains(int id) {
        for (int i = 0; i < size(); i++) {
            if (get(i) == id) {
                return(true);
            }
        }
        return(false);
    }


    void sortvec() {
        sort(m_index.begin(), m_index.end());
    }


private:
    svec<int> m_index;
    int m_id;

};


#endif

