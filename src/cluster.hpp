
#ifndef CLUSTER_HPP
#define CLUSTER_HPP

#include "pmdf.hpp"
#include "ovl.hpp"
//#include <random>
#include <iostream>
#include <set>

template <typename T, class P> class Cluster
{
    private:
        int clusterId;
        vector<T> centroid_pos; // the positions of centroid, for plotting
        vector<P> pmdfs;

        vector<T> setUnion(vector<T>& v1, vector<T>& v2){

          set<T> s1(v1.begin(), v1.end());
          set<T> s2(v2.begin(), v2.end());
          set<T> intersect;
          set_union(s1.begin(), s1.end(), s2.begin(), s2.end(),
                           std::inserter(intersect, intersect.begin()));
          std::vector<T> v(intersect.begin(), intersect.end());
          return v;

        }

    public:

        Cluster(){}

        Cluster(int clusterId)  // Partial construction, be used in kmeans.
        {
            this->clusterId = clusterId;
        }

        Cluster(int clusterId, vector<P> const &pmdfs) // full construction, for testing.
        {
            this->clusterId = clusterId;
            this->pmdfs = pmdfs;
        }

        ~Cluster(){}

        void removeAllPmdfs() { pmdfs.clear(); }

        int getId() { return clusterId; }

        void setId(int id) {this->clusterId = id;}

        P &getPmdf(int pos) { return pmdfs[pos]; }

        int getSize() { return pmdfs.size(); }

        vector<int> getPmdfIds()
        {
            vector<int> pmdfIds;
            for (int i = 0; i < getSize(); ++i)
            {
                pmdfIds.push_back(pmdfs[i].getID());
            }
            return pmdfIds;
        }

        vector<P> &getPmdfs() {return this->pmdfs;}

        void setPmdfs(vector<P> &pmdfs) {this->pmdfs = pmdfs;}

        void addPmdf(P &p)
        {
            p.setCluster(this->clusterId);
            pmdfs.push_back(p);
        }

        bool removePmdf(int pId)
        {
            int size = pmdfs.size();

            for (int i = 0; i < size; i++)
            {
                if (pmdfs[i].getID() == pId)
                {
                    pmdfs.erase(pmdfs.begin() + i);
                    return true;
                }
            }
            return false;
        }

        T getMinVal()
        {
          int n = this->getSize();
          T min_val = pmdfs[0].getMinVal();

          for (int i = 1; i < n; ++i)
          {
            if(min_val>pmdfs[i].getMinVal()){
              min_val = pmdfs[i].getMinVal();
            }
          }
          return min_val;
        }

        T getMaxVal()
        {
          int n = this->getSize();
          T max_val = pmdfs[0].getMaxVal();

          for (int i = 1; i < n; ++i)
          {
            if(max_val<pmdfs[i].getMaxVal()){
              max_val = pmdfs[i].getMaxVal();
            }
          }
          return max_val;
        }
        // initialize centroid positions
        void initCentroidPos(int grid)
        {
          bool type = pmdfs[0].getType();
          if(type){
            T min_val = getMinVal();
            T max_val = getMaxVal();
            double step = double((max_val-min_val)/(grid-1));
            for (int i = 0; i < grid; ++i)
            {
              this->centroid_pos.push_back(min_val+i*step);
            }
          }else{
            int n = this->getSize();
            vector<T> vec = pmdfs[0].getVals();
            for (int i = 0; i < n; ++i){
              vec = setUnion(vec,pmdfs[i].getVals());
            }
            this->centroid_pos = vec;
          }
        }
        // get centroid positions
        vector<T> &getCentroidPos()
        {
          return this->centroid_pos;
        }

        // for plot, return the probability at position pos.
        double centroidEstimator(T pos)
        {
          bool type = pmdfs[0].getType();
          double tpres=0;
          int n = this->getSize();

          for (int i = 0; i < n; ++i)
          {
            if(type) // is pdf
            {
              kdepp::Kde1d<T> densEst = pmdfs[i].getDensityEstimator();
              tpres = tpres + densEst.eval(pos);
            }else{  // is pmf
              // there are n pmfs
              int count = 0;
              int pmf_size = pmdfs[i].getSize();
              for (int j = 0; j < pmf_size; ++j)
              {
                if(pos==pmdfs[i].getVal(j)){ // count number pos in pmf
                  count = count + 1;
                }
              }
              tpres = tpres + double(count)/double(pmf_size);
            }
          }
          return tpres/n;
        }
        // calculate centroid at positions
        vector<double> calculateCentroid(int grid=1000)
        {
          vector<double> results;
          initCentroidPos(grid);
          int n=centroid_pos.size();
          for (int i = 0; i < n; ++i)
          {
            results.push_back(centroidEstimator(centroid_pos[i]));
          }
          return results;
        }

        // calculate Am
        double calculateAm(int grid=1000){
          double Am = 0;
          OVL ovl;
          for (int i = 0; i < this->getSize(); i++){
            Am = Am + ovl.ovlEst(pmdfs, pmdfs[i], grid);
          }
          return Am;
        }

};

// // calculate g_m(x).
// void calculateCentroid(double bootstrapScale)
// {
//
//     bool type = pmdfs[0].getType();  // 1: double, pdf; 0: int, pmf
//
//     int cSize = this->getSize(); // n_m
//
//     int sumSize = 0;
//
//     for (int i = 0; i < cSize; ++i)
//     {
//           sumSize = sumSize+pmdfs[i].getSize();
//     }
//     int bootstrapN;
//     if(sumSize<1000){
//       bootstrapN=1000;
//     }else{
//       if(cSize==1){
//         bootstrapN = int(sumSize*bootstrapScale);
//       }else if(cSize==2){
//         bootstrapN = int(sumSize*bootstrapScale*0.5);
//       }else{
//         bootstrapN = sumSize;
//       }
//     }
//
//     vector<T> c_pmdf_vals; // centroid values
//
//     if(type){// continuous distribution
//       // 1. create positions
//       // (1) find range
//       T min_val, max_val;
//       max_val = pmdfs[0].getMaxVal();
//       min_val = pmdfs[0].getMinVal();
//       for(int i = 1; i < cSize; ++i){
//         if(min_val>pmdfs[i].getMinVal()){
//           min_val = pmdfs[i].getMinVal();
//         }
//         if(max_val<pmdfs[i].getMaxVal()){
//           max_val = pmdfs[i].getMaxVal();
//         }
//       }
//       // (2) create positions in above range
//       bootstrapN=50;
//       double step = double((max_val-min_val)/(bootstrapN-1));
//       // (3)
//       for (int i = 0; i < bootstrapN; ++i)
//       {
//         double temp = 0;
//         for(int j = 0; j < cSize; ++j){
//           kdepp::Kde1d<T> kdestor = pmdfs[j].getDensityEstimator();
//           temp = temp + kdestor.eval(min_val+i*step);
//           cout<<temp<<endl;
//         }
//         c_pmdf_vals.push_back(temp/cSize);
//       }
//
//     }else{ // discrete distribution
//       std::random_device rd;  //Will be used to obtain a seed for the random number engine
//       std::mt19937 gen(rd()); //Standard mersenne_twister_engine seeded with rd()
//       std::uniform_int_distribution<> distUnif(0, cSize-1); // the number of pmdfs
//       for (int i = 0; i < bootstrapN; ++i)
//       {
//         int pmdf_id = distUnif(gen);
//         int rd = rand() % pmdfs[pmdf_id].getSize();
//         c_pmdf_vals.push_back(pmdfs[pmdf_id].getVals()[rd]);
//       }
//     }
//     centroid.setVals(c_pmdf_vals);
// }

#endif
