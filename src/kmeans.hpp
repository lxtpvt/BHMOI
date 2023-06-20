
#ifndef KMEANS_HPP
#define KMEANS_HPP

#include <omp.h>
#include "pmdf.hpp"
#include "cluster.hpp"
#include "ovl.hpp"
#include <Rcpp.h>
#include <cmath>
#include<cstdlib>

using namespace Rcpp;

template <class C, class P> class KMeans
{
    private:
        int K, iters, total_pmdfs;
        vector<C> clusters;

        void clearClusters()
        {
            for (int i = 0; i < K; i++)
            {
                clusters[i].removeAllPmdfs();
            }
        }

        double getDistance2Cluster(C &cluster, P &pmdf, int grid=1000)
        {
          OVL ovl;
          return ovl.ovlDist(cluster.getPmdfs(), pmdf, grid);
        }

        int getNearestClusterId(P &pmdf, int n, bool weighted, double b, int grid=1000)
        {
            double min_dist;
            int NearestClusterId;
            if(weighted){
              min_dist = getDistance2Cluster(clusters[0], pmdf, grid)*pow((1-clusters[0].getSize()/n), b);
            }else{
              min_dist = getDistance2Cluster(clusters[0], pmdf, grid);
            }
            NearestClusterId = clusters[0].getId();

            for (int i = 1; i < K; i++)
            {
              //Rcout << "middle, K : " << i << std::endl;
                double dist;
                if(weighted){
                  dist = getDistance2Cluster(clusters[i], pmdf, grid)*pow((1-clusters[i].getSize()/n), b);
                }else{
                  dist = getDistance2Cluster(clusters[i], pmdf, grid);
                }
                if (dist < min_dist)
                {
                    min_dist = dist;
                    NearestClusterId = clusters[i].getId();
                }
            }
            return NearestClusterId;
        }

    public:
        KMeans(int K, int iterations)
        {
            this->K = K;
            this->iters = iterations;
        }

        KMeans(){}

        ~KMeans(){}

        int getK() {return this->K;}

        void setK(int k) {this->K = k;}

        int getIters() {return this->iters;}

        void setIters(int iters) {this->iters = iters;}

        int getSizePmdfs() {return this->total_pmdfs;}

        vector<C> &getClusters() {return this->clusters;}

        void initClusters(vector<P> &pmdfs, int seed) // Initializing Clusters
        {
            total_pmdfs = pmdfs.size();

            vector<int> used_pIds;

            // K is the number of clusters. Randomly choose K pmdfs as the centroid of K clusters.
            for (int i = 1; i <= K; i++) //clusterIDs are started from 1.
            {
                while (true)
                {
                    srand(seed+i);
                    //srand(seed);
                    int index = rand() % total_pmdfs;

                    if (find(used_pIds.begin(), used_pIds.end(), index) ==
                        used_pIds.end())
                    {
                        used_pIds.push_back(index);
                        pmdfs[index].setCluster(i);
                        C cluster(i); // using the partial construction
                        cluster.addPmdf(pmdfs[index]);
                        clusters.push_back(cluster);
                        break;
                    }
                }
            }
            //Rcout << "Clusters initialized = " << clusters.size() << std::endl;
        }

        void run(vector<P> &pmdfs,bool weighted, double b, int nthreads, int seed, int grid)  // run kmeans
        {
            total_pmdfs = pmdfs.size();

            initClusters(pmdfs,seed); // Initializing Clusters

            int iter = 1;
            while (true)
            {
                bool done = true;

                #pragma omp parallel for num_threads(nthreads)
                // Add all pmdfs to their nearest cluster. for each pmdf, find its nearest cluster.
                for (int i = 0; i < total_pmdfs; i++)
                {
                    int currentClusterId = pmdfs[i].getCluster();
                    int nearestClusterId = getNearestClusterId(pmdfs[i],total_pmdfs, weighted, b, grid);
                    if (currentClusterId != nearestClusterId)
                    {
                        pmdfs[i].setCluster(nearestClusterId);
                        done = false;
                    }
                } // all pmdfs be assigned clusterId.

                clearClusters(); // clear all existing clusters

                // all clusters be reassigned pmdfs
                for (int i = 0; i < total_pmdfs; i++)
                {
                  clusters[pmdfs[i].getCluster() - 1].addPmdf(pmdfs[i]);
                }

                if (done || iter >= iters)
                {
                    break;
                }
                iter++;
            } // end of while loop. all the clustering information can be found in class member clusters

        } // end or run
};

#endif
