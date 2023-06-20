
/* The coefficient of overlap for any pair of pdf or pmf*/

#ifndef OVL_HPP
#define OVL_HPP

#include "pmdf.hpp"
#include "kde.hpp" // For estimating empirical density.
#include "cluster.hpp"
#include <algorithm> //std::sort
#include <set>

class OVL
{
    private:

    	double scaleBW;
    	int sizeT;

    	vector<int> intersection(vector<int>& v1, vector<int>& v2){
    	  set<int> s1(v1.begin(), v1.end());
    	  set<int> s2(v2.begin(), v2.end());
    	  set<int> intersect;
    	  set_intersection(s1.begin(), s1.end(), s2.begin(), s2.end(),
                        std::inserter(intersect, intersect.begin()));
    	  std::vector<int> v(intersect.begin(), intersect.end());
    	  return v;
    	}

	public:

        OVL(double bw, int size)
        {
        	this->scaleBW = bw;
        	this->sizeT = size;
        }

        OVL()
        {
        	this->scaleBW = 0.01;
        	this->sizeT = 75;
        }

        ~OVL(){}

        double getScaleBW(){return this->scaleBW;}

        void setScaleBW(double scale){this->scaleBW = scale;}

        int getSizeT(){return this->sizeT;}

        void setSizeT(int size){this->sizeT = size;}

	      //=============================
		    // if they are two pdfs
        double ovlEst(Pmdf<double> &pdf1, Pmdf<double> &pdf2, int grid = 1000)
        {
        	int n,m,smn;
        	double l, r, rg, lvals, step;

        	n = pdf1.getSize();
        	m = pdf2.getSize();
        	l = min(pdf1.getMinVal(),pdf2.getMinVal());
        	r = max(pdf1.getMaxVal(),pdf2.getMaxVal());
        	rg = r - l;
        	lvals = l - 0.1*rg;
        	//rvals = r + 0.2*rg;
        	step = rg*1.2/(grid-1);

        	// create a vector values with equal step
    			std::vector<double> dvec(grid);
    			for (int i = 0; i < grid; ++i)
    			{
    				dvec[i] = lvals + i*step;
    			}

        	smn = min(n,m);
        	kdepp::Kde1d<double> densEst1 = pdf1.getDensityEstimator();
        	kdepp::Kde1d<double> densEst2 = pdf2.getDensityEstimator();
        	double bw1 = densEst1.get_bandwidth();
        	double bw2 = densEst2.get_bandwidth();
        	densEst1.set_bandwidth(bw1*scaleBW);
        	densEst2.set_bandwidth(bw2*scaleBW);

        	if (smn<75) {
        		// do delta-1
        		double dSum = 0;
        		double dSum1 = 0;
        		double dSum2 = 0;
        		vector<double> densVec1,densVec2;

        		for (int i = 0; i < grid; ++i)
        		{
        		  double tp1 = densEst1.eval(dvec[i]);
        		  double tp2 = densEst2.eval(dvec[i]);
        		  densVec1.push_back(tp1);
        		  densVec2.push_back(tp2);
        			dSum1 = dSum1 + tp1;
        			dSum2 = dSum2 + tp2;
        		}

        		for (int i = 0; i < grid; ++i)
        		{
        			dSum = dSum + min(densVec1[i]/dSum1, densVec2[i]/dSum2);
        		}

        		return dSum;

        	} else {
        		// do detal-4
        		double dSum1 = 0;
        		for (int i = 0; i < n; ++i)
        		{
        			dSum1 = dSum1 + min(1.0,densEst2.eval(pdf1.getVal(i))/densEst1.eval(pdf1.getVal(i)));
        		}
        		double dSum2 = 0;
        		for (int i = 0; i < m; ++i)
        		{
        			dSum2 = dSum2 + min(1.0,densEst1.eval(pdf2.getVal(i))/densEst2.eval(pdf2.getVal(i)));
        		}
        		return 0.5*(dSum1/n + dSum2/m);

        	}
        }
	      //=============================
	      // if they are two pmfs
    	  double ovlEst(Pmdf<int> &pmf1, Pmdf<int> &pmf2)
    	  {
    	    vector<int> d1 = pmf1.getVals();
    	    vector<int> d2 = pmf2.getVals();
    	    vector<int> intersect_e = intersection(d1, d2);
    	    double res=0;
    	    for (vector<int>::iterator e = intersect_e.begin(); e != intersect_e.end(); ++e){
    	      int count1=0;
    	      int count2=0;
    	      for (vector<int>::iterator e1 = d1.begin(); e1 != d1.end(); ++e1){
    	        if(*e==*e1){
    	          count1 = count1 + 1;
    	        }
    	      }
    	      for (vector<int>::iterator e2 = d2.begin(); e2 != d2.end(); ++e2){
    	        if(*e==*e2){
    	          count2 = count2 + 1;
    	        }
    	      }
    	      double p1 = double(count1)/d1.size();
    	      double p2 = double(count2)/d2.size();
    	      if(p1<p2){
    	        res=res+p1;
    	      }else{
    	        res=res+p2;
    	      }
    	    }
    	    return res;
    	  }
	      //=============================
	      // if they are a cluster and a pdf or pmf
	      double ovlEst(vector<Pmdf<double>> &cluster, Pmdf<double> &pmdf, int grid=1000) // it is pdf
	      {
	        int nPmdfs = cluster.size();
	        double sumOVL=0;
	        for (int i = 0; i < nPmdfs; ++i)
	        {
	          sumOVL = sumOVL + ovlEst(cluster[i],pmdf,grid);
	        }
	        return sumOVL/nPmdfs;
	      }

    	  double ovlEst(vector<Pmdf<int>> &cluster, Pmdf<int> &pmdf, int grid=1000)  // it is pmf
    	  {
    	    int nPmdfs = cluster.size();
    	    double sumOVL=0;
    	    for (int i = 0; i < nPmdfs; ++i)
    	    {
    	      sumOVL = sumOVL + ovlEst(cluster[i],pmdf);
    	    }
    	    return sumOVL/nPmdfs;
    	  }
	      //=======================================================================
        // Define distance by using OVL
        //pdf, need grid
        double ovlDist(Pmdf<double> &pmdf1, Pmdf<double> &pmdf2, int grid = 1000) // it is pdf
        {
        	return (1.0-ovlEst(pmdf1,pmdf2,grid));
        }

    	  double ovlDist(Pmdf<int> &pmdf1, Pmdf<int> &pmdf2, int grid = 1000) // it is pmf
    	  {
    	    return (1.0-ovlEst(pmdf1,pmdf2));
    	  }

    	  double ovlDist(vector<Pmdf<double>> &cluster, Pmdf<double> &pmdf, int grid = 1000) // it is pdf
    	  {
    	    return (1.0-ovlEst(cluster,pmdf,grid));
    	  }

    	  double ovlDist(vector<Pmdf<int>> &cluster, Pmdf<int> &pmdf, int grid = 1000) // it is pmf
    	  {
    	    return (1.0-ovlEst(cluster,pmdf));
    	  }
};

#endif
