
#ifndef PMDF_HPP
#define PMDF_HPP

//#include <iostream>
#include <vector>
#include "kde.hpp" // For estimating empirical density.

using namespace std;

template <typename T> class Pmdf
{
    private:

        int pId, clusterId, size;
        vector<T> vals;
        bool type; // 1: pdf; 0: pmf

    public:

        Pmdf(int id, vector<T> const &v, bool t)
        {
            this->pId = id;
            this->size = v.size();
            this->vals = v;
            // Initially not assigned to any cluster
            this->clusterId = -1;
            this->type = t;
        }

        Pmdf(){}

        ~Pmdf(){}
      
        bool getType() { return this->type; }

        int getID() { return pId; }

        void setID(int id) {this->pId = id;}

        int getSize() { return size; }

        int getCluster() { return clusterId; }

        void setCluster(int cId) { clusterId = cId; }

        T getVal(int pos) { return vals[pos]; }

        void setVal(int pos, T val) {vals[pos] = val;}

        vector<T> &getVals(){return vals;}

        void setVals(vector<T> &vals)
        {
            this->vals = vals;
            this->size = vals.size();
        }

        T getMaxVal(){return *max_element(vals.begin(), vals.end()); }

        T getMinVal(){return *min_element(vals.begin(), vals.end()); }

        // density estimator
        kdepp::Kde1d<T> getDensityEstimator()
        {
            kdepp::Kde1d<T> kernel(this->vals);
            return kernel;
        }

        kdepp::Kde1d<T> getDensityEstimator(double h)
        {
            kdepp::Kde1d<T> kernel(this->vals);
            kernel.set_bandwidth(h);
            return kernel;
        }

};

#endif