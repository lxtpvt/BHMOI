#include "kmeans.hpp"
#include "pmdf.hpp"
#include "cluster.hpp"
#include "ovl.hpp"
#include "kde.hpp"
#include "progress.hpp"
#include "progress_bar.hpp"
#include <Rcpp.h>

using namespace Rcpp;

// [[Rcpp::export(.kde)]]
NumericVector docbi_kde(
	NumericVector data,			// data
	NumericVector positions	// estimating positions
)
{
	std::vector<double> vals = as<std::vector<double> >(data);

	Pmdf<double> pdf(1,vals,1);

	kdepp::Kde1d<double> kdestor = pdf.getDensityEstimator();
	//kdestor.set_bandwidth(kdestor.get_bandwidth());

	NumericVector res;
	for (int i = 0; i < positions.length(); ++i)
	{
	  double temp = kdestor.eval(positions[i]);
		res.push_back(temp);
	}
	return res;
}

// [[Rcpp::export(.clusterKDE)]]
NumericVector docbi_clusterKDE(
    List data,			// data
    NumericVector positions	// estimating positions
)
{
  Cluster<double,Pmdf<double>> cluster(1);
  for (int i = 0; i < data.length(); ++i)
  {
    std::vector<double> vals = as<std::vector<double> >(data[i]);
    Pmdf<double> pdf(i,vals,1);
    cluster.addPmdf(pdf);
  }

  NumericVector res;
  for (int i = 0; i < positions.length(); ++i)
  {
    double temp = cluster.centroidEstimator(positions[i]);
    res.push_back(temp);
  }
  return res;
}


// [[Rcpp::export(.ovl_continuous)]]
double docbi_ovl_continuous( NumericVector data_1, NumericVector data_2, size_t grid)
{
  std::vector<double> d1 = as<std::vector<double> >(data_1);
  std::vector<double> d2 = as<std::vector<double> >(data_2);
  Pmdf<double> pdf1(1,d1,1);
  Pmdf<double> pdf2(2,d2,1);

  OVL ovlEstor;
  return ovlEstor.ovlEst(pdf1,pdf2,grid);
}

// [[Rcpp::export(.ovl_discrete)]]
double docbi_ovl_discrete( NumericVector data_1, NumericVector data_2)
{
  std::vector<int> d1 = as<std::vector<int> >(data_1);
  std::vector<int> d2 = as<std::vector<int> >(data_2);
  Pmdf<int> pmf1(1,d1,0);
  Pmdf<int> pmf2(2,d2,0);

  OVL ovlEstor;
  return ovlEstor.ovlEst(pmf1,pmf2);
}

// [[Rcpp::export(.KMeanCluster)]]
List docbi_kmeansCluster( size_t K, List dataSets, bool weighted, double b, size_t n_sim, size_t nthreads,
                          size_t nIters, bool showmessage, bool is_pdf, int seed, int grid=1000)
{
  if(showmessage){
    Rcout << "Running Weighted K-Means Clustering..." << std::endl;
    Rcout << "The number of clusters, K = " << K << std::endl;
  }
  List resList;
  Progress p(n_sim, showmessage);
  for (int j = 0; j < n_sim; ++j)
  {
    List Clusters;
    List Centroids;
    List Ams;
    if(is_pdf){ // if it's pdf
      KMeans<Cluster<double,Pmdf<double>>, Pmdf<double>> kmc_pdf(K, nIters);
      vector<Pmdf<double>> pdfs;
      for (int i = 0; i < dataSets.length(); ++i)
      {
        std::vector<double> tpv = as<std::vector<double> >(dataSets[i]);
        Pmdf<double> temp_pmdf(i,tpv,1);
        pdfs.push_back(temp_pmdf);
      }
      kmc_pdf.run(pdfs,weighted,b,nthreads,seed,grid);
      for (int i = 0; i < K; ++i)
      {
        vector<int> ids = kmc_pdf.getClusters()[i].getPmdfIds();
        transform(ids.begin(), ids.end(), ids.begin(),
                  bind2nd(std::plus<int>(), 1));
        Clusters.push_back(ids);
        List tpCentroid;
        vector<double>cprobs = kmc_pdf.getClusters()[i].calculateCentroid();
        vector<double>cpos = kmc_pdf.getClusters()[i].getCentroidPos();
        tpCentroid.push_back(cpos);
        tpCentroid.push_back(cprobs);
        Centroids.push_back(tpCentroid);
        Ams.push_back(kmc_pdf.getClusters()[i].calculateAm(grid));
      }
    }else{
      KMeans<Cluster<int, Pmdf<int>>, Pmdf<int>> kmc_pmf(K, nIters);
      vector<Pmdf<int>> pmfs;
      for (int i = 0; i < dataSets.length(); ++i)
      {
        std::vector<int> tpv = as<std::vector<int> >(dataSets[i]);
        Pmdf<int> temp_pmdf(i,tpv,0);
        pmfs.push_back(temp_pmdf);
      }
      kmc_pmf.run(pmfs,weighted,b,nthreads,seed,1);
      for (int i = 0; i < K; ++i)
      {
        vector<int> ids = kmc_pmf.getClusters()[i].getPmdfIds();
        transform(ids.begin(), ids.end(), ids.begin(),
                  bind2nd(std::plus<int>(), 1));
        Clusters.push_back(ids);
        List tpCentroid;
        vector<double>cprobs = kmc_pmf.getClusters()[i].calculateCentroid();
        vector<int>cpos = kmc_pmf.getClusters()[i].getCentroidPos();
        tpCentroid.push_back(cpos);
        tpCentroid.push_back(cprobs);
        Centroids.push_back(tpCentroid);
        Ams.push_back(kmc_pmf.getClusters()[i].calculateAm(grid));
      }
    }
    resList.push_back(Rcpp::List::create(Clusters,Ams,Centroids));
    p.increment();
  }

  return resList;
}
