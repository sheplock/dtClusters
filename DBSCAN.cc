#include "DBSCAN.h"
#include "TMath.h"
#include <iostream>
#include "TVector3.h"
#include "TGraph.h"
#include "TF1.h"

struct largest_nDt_cluster_
{
  inline bool operator() (const dtCluster& c1, const dtCluster& c2){return c1.nDtSegments > c2.nDtSegments;}
} largest_nDt_cluster;

struct hits
{
  float time;
  float error;
  bool strip;
};
const double theWireError_ = 8.6;
const double theStripError_ = 7.0;
const double thePruneCut_ = 9.0;

// function called to run the DBSCAN algorithm
int DBSCAN::run()
{
    int clusterID = 1; // index used to classify which cluster each point belongs to
    vector<Point>::iterator iter; //iterator for DT hits
    for(iter = m_points.begin(); iter != m_points.end(); ++iter)
    {
        if ( iter->clusterID == UNCLASSIFIED ) // check point hasn't been classified yet in another iteration of algorithm
        {
            // result is failure if point can't seed new cluster
            // if point isn't failure (does seed a new cluster) then increase clusterID to start new cluster
            if ( expandCluster(*iter, clusterID) != FAILURE )
            {
                clusterID += 1;
            }
        }
    }
    // total number of clusters is clusterID-1 (initialized to 1 before any clusters made)
    nClusters = clusterID-1;
    return (clusterID-1);
}

int DBSCAN::expandCluster(Point point, int clusterID)
{
    //calculate all neighbors within distance parameter, including point itself
    vector<int> clusterSeeds = calculateCluster(point);

    // if fewer neighboring points than min cluster size, point is labeled as noise
    if ( clusterSeeds.size() < m_minPoints )
    {
        point.clusterID = NOISE;
        return FAILURE;
    }
    else // point is a core point
    {
        int index = 0, indexCorePoint = 0;
        vector<int>::iterator iterSeeds;
        // loop through neighbors of core points
        for( iterSeeds = clusterSeeds.begin(); iterSeeds != clusterSeeds.end(); ++iterSeeds)
        {
            // setting all the neighbors to the same clusterID
            m_points.at(*iterSeeds).clusterID = clusterID;
            // get the index of the core point itself
            if (m_points.at(*iterSeeds).eta == point.eta && m_points.at(*iterSeeds).phi == point.phi )
            {
                indexCorePoint = index;
            }
            // increment index for each time iterating through the loop
            ++index;
        }
        // remove the core point
        clusterSeeds.erase(clusterSeeds.begin()+indexCorePoint);
        // now clusterSeeds only contains neighbors, not core point itself
        // loop through the neighbors again
        // if one of the neighbors is itself a core point, all points in the epsilon neighborhood of that point will be added to cluster
        for( vector<int>::size_type i = 0, n = clusterSeeds.size(); i < n; ++i )
        {
            // check if the neighbors themselves are core points
            vector<int> clusterNeighors = calculateCluster(m_points.at(clusterSeeds[i]));
            if ( clusterNeighors.size() >= m_minPoints ) // this neighbor point is a core point
            {
                vector<int>::iterator iterNeighors;
                // iterate through the neighbors of the new core point
                for ( iterNeighors = clusterNeighors.begin(); iterNeighors != clusterNeighors.end(); ++iterNeighors )
                {
                    // new points that don't already belong to a cluster
                    if ( m_points.at(*iterNeighors).clusterID == UNCLASSIFIED || m_points.at(*iterNeighors).clusterID == NOISE )
                    {
                        // if point hasn't been classified yet, add it to clusterSeeds so it will be checked for being a core point
                        // if noise, already know point isn't a core point from a previous iteration
                        // if point is added, change n so the loop through all clusterSeeds contains new point
                        if ( m_points.at(*iterNeighors).clusterID == UNCLASSIFIED )
                        {
                            clusterSeeds.push_back(*iterNeighors);
                            n = clusterSeeds.size();
                        }
                        // label points as belonging to this cluster
                        m_points.at(*iterNeighors).clusterID = clusterID;
                    }
                }
            }
        }

        return SUCCESS;
    }
}

// get index of all points in epsilon neighborhood of point
vector<int> DBSCAN::calculateCluster(Point point)
{
    int index = 0;
    vector<Point>::iterator iter;
    vector<int> clusterIndex;
    Point minimum;
    // iterate through all points, store index of points within epsilon
    for( iter = m_points.begin(); iter != m_points.end(); ++iter)
    {
        if ( calculateDistance(point, *iter) <= m_epsilon )
        {
            clusterIndex.push_back(index);
        }
        index++;
    }
    return clusterIndex;
}

inline double DBSCAN::calculateDistance( Point pointCore, Point pointTarget )
{
    // return sqrt(pow(pointCore.x - pointTarget.x,2)+pow(pointCore.y - pointTarget.y,2)+pow(pointCore.z - pointTarget.z,2));
    // return sqrt(pow(pointCore.eta - pointTarget.eta,2)+pow(deltaPhi(pointCore.phi, pointTarget.phi),2)+pow((pointCore.t - pointTarget.t)/100.0,2));
    return sqrt(pow(pointCore.eta - pointTarget.eta,2)+pow(deltaPhi(pointCore.phi, pointTarget.phi),2));

}

void DBSCAN::clear_clusters(){
  // nClusters = 0;
  clusterSize.clear();
  dtLabels.clear();
  clusterEta.clear();
  clusterPhi.clear();
  clusterX.clear();
  clusterY.clear();
  clusterZ.clear();
  clusterTime.clear();
  clusterTimeTotal.clear();
  clusterTimeTotal.clear();
  clusterMajorAxis.clear();
  clusterMinorAxis.clear();
  clusterXSpread.clear();
  clusterYSpread.clear();
  clusterXYSpread.clear();
  clusterRSpread.clear();
  clusterZSpread.clear();
  clusterTimeSpread.clear();
  clusterTimeSpreadWeighted.clear();
  clusterTimeSpreadWeightedAll.clear();

  clusterEtaPhiSpread.clear();
  clusterEtaSpread.clear();
  clusterPhiSpread.clear();
  clusterDeltaRSpread.clear();

  dtCluster.clear();
}
int DBSCAN::result(){

  // for (unsigned int i = 0;i < m_pointSize;i++)
  // {
  //   dtLabels.push_back(m_points[i].clusterID);
  //
  // }
  for(int i = 0; i < nClusters; i++)
  {
    float avg_x(0.0), avg_y(0.0), avg_z(0.0), avg_tWire(0.0), avg_tWirePruned(0.0), avg_t(0.0), avg_tTotal(0.0),tTotalSpreadPruned(0.0);
    float avg_x_sl2(0.0), avg_y_sl2(0.0), avg_z_sl2(0.0);
    float avg_eta(0.0), avg_phi(0.0);
    int size(0), size_z(0), size_xy(0);
    vector<float> wireTimes;
    vector<float> stripTimes;
    std::vector<hits> dtHits;

    vector<Point>::iterator iter;
    for(iter = m_points.begin(); iter != m_points.end(); ++iter)
    {
      if ( iter->clusterID == i+1 )
      {
          if (iter->superlayer == 2) //for DT rechits that only have coordinates in Z
          {
            avg_x_sl2 += iter->x;
            avg_y_sl2 += iter->y;
            avg_z_sl2 += iter->z;
            size_z++;
          }
          else if (iter->superlayer == 1 || iter->superlayer == 3)
          {
            avg_x += iter->x;
            avg_y += iter->y;
            avg_z += iter->z;
            avg_t += iter->t;
            size_xy ++;
          }
          else //csc or for DT "wrong" rechit coordinates
          {
            avg_x += iter->x;
            avg_y += iter->y;
            avg_z += iter->z;
            avg_t += iter->t;
            avg_tWire += iter->twire;
            wireTimes.push_back(iter->twire);
            stripTimes.push_back(iter->t);
            hits thisHit;
            thisHit.time = iter->twire;
            thisHit.error = 1./(theWireError_*theWireError_);
	          thisHit.strip = false;
            dtHits.push_back(thisHit);
            thisHit.time = iter->t;
            thisHit.error = 1./(theStripError_*theStripError_);
	          thisHit.strip = true;
            dtHits.push_back(thisHit);

          }
          size ++;

      }
    }
    if (size_xy > 0 && size_z > 0) //for DT correct position, calculate average Z using sl2 and average XY using sl1/3
    {
      avg_x = avg_x/size_xy;
      avg_y = avg_y/size_xy;
      avg_z = avg_z_sl2/size_z;
    }
    else if (size_xy == 0 && size_z == 0) //csc or DT wrong position
    {
      avg_x = avg_x/size;
      avg_y = avg_y/size;
      avg_z = avg_z/size;
      // cout<<avg_x<<","<<avg_y<<","<<avg_z<<endl;
    }
    else if (size_xy > 0 && size_z == 0)
    {
      avg_x = avg_x/size_xy;
      avg_y = avg_y/size_xy;
      avg_z = avg_z/size_xy;

    }
    else
    {
      avg_x = avg_x_sl2/size_z;
      avg_y = avg_y_sl2/size_z;
      avg_z = avg_z_sl2/size_z;

    }
    avg_tTotal = (avg_t + avg_tWire)/(2 * size);
    avg_t = avg_t/size;
    avg_tWire = avg_tWire/size;

    // prune wire time
    //The wire times have a long tail that has to be pruned.  The strip times (tpeak) are fine
    bool modified = true;
    while (modified) {
      modified = false;
      double maxDiff = -1;
      std::vector<float>::iterator maxHit;
      for (std::vector<float>::iterator itWT = wireTimes.begin(); itWT != wireTimes.end(); ++itWT) {
        float diff = fabs(*itWT - avg_tTotal);
        if (diff > maxDiff) {
          maxDiff = diff;
          maxHit = itWT;
        }
      }
      if (maxDiff > 26) {
        int N = size + wireTimes.size();
        avg_tTotal = (avg_tTotal * N - (*maxHit)) / (N - 1);
        wireTimes.erase(maxHit);
        modified = true;
      }
    }

    //new timing calculation, error weighted
    // https://github.com/cms-sw/cmssw/blob/master/RecoMuon/MuonIdentification/src/CSCTimingExtractor.cc
    modified = false;
    double totalWeightTimeVtx = 0;
    double timeVtx = 0;
    double timeSpread = 0;
    do {
      modified = false;
      totalWeightTimeVtx = 0;
      timeVtx = 0;
      timeSpread = 0;
      for (std::vector<hits>::iterator it = dtHits.begin(); it != dtHits.end(); ++it) {
        timeVtx += it->time * it->error;
        totalWeightTimeVtx += it->error;
      }
      timeVtx /= totalWeightTimeVtx;

      // cut away outliers
      double diff_tvtx;
      double chimax = 0.0;
      int tmmax;
      for (unsigned int i = 0; i < dtHits.size(); i++) {
        diff_tvtx = (dtHits[i].time - timeVtx) * (dtHits[i].time - timeVtx) * dtHits[i].error;

        if (diff_tvtx > chimax) {
          tmmax =  i;
          chimax = diff_tvtx;
        }
      }
      // cut away the outliers
      if (chimax > thePruneCut_) {
        dtHits.erase(dtHits.begin()+tmmax);
        modified = true;
      }
    } while (modified);
    int count = 0;
    for (std::vector<hits>::iterator it = dtHits.begin(); it != dtHits.end(); ++it) {
      if (it->strip)
      {
        timeSpread += (it->time - timeVtx)*(it->time - timeVtx);
        count++;
      }

    }
    timeSpread = sqrt(timeSpread/count);


    // calculate cluster eta and phi
    avg_phi = atan(avg_y/avg_x);
    if  (avg_x < 0.0){
      avg_phi = TMath::Pi() + avg_phi;
    }
    avg_phi = deltaPhi(avg_phi,0.0);
    avg_eta = atan(sqrt(pow(avg_x,2)+pow(avg_y,2))/abs(avg_z));
    avg_eta = -1.0*TMath::Sign(1.0, avg_z)*log(tan(avg_eta/2));

    clusterEta.push_back(avg_eta);
    clusterPhi.push_back(avg_phi);
    clusterX.push_back(avg_x);
    clusterY.push_back(avg_y);
    clusterZ.push_back(avg_z);
    clusterTime.push_back(avg_t);
    clusterTimeTotal.push_back(avg_tTotal);
    clusterSize.push_back(size);

    clusterTimeWeighted.push_back(timeVtx);
    clusterTimeSpreadWeighted.push_back(timeSpread);

  }
  return 0;
}
int DBSCAN::clusterMoments()
{

  for(int i = 0; i < nClusters; i++)
  {
    float m11(0.0), m12(0.0), m22(0.0);
    float XSpread(0.0), YSpread(0.0), ZSpread(0.0), TSpread(0.0),  TSpreadAll(0.0), XYSpread(0.0), RSpread(0.0), DeltaRSpread(0.0);



    vector<Point>::iterator iter;
    for(iter = m_points.begin(); iter != m_points.end(); ++iter)
    {
      if ( iter->clusterID == i+1 )
      {

          m11 += (iter->eta-clusterEta[i])*(iter->eta-clusterEta[i]);
          m12 += (iter->eta-clusterEta[i])* deltaPhi(iter->phi,clusterPhi[i]);
          m22 += deltaPhi(iter->phi,clusterPhi[i])*deltaPhi(iter->phi,clusterPhi[i]);
          DeltaRSpread +=  pow(deltaR(clusterEta[i], clusterPhi[i], iter->eta, iter->phi),2);
          XYSpread += (iter->x - clusterX[i])*(iter->y - clusterY[i]);
          XSpread += (iter->x - clusterX[i]) * (iter->x - clusterX[i]);
          YSpread += (iter->y - clusterY[i]) * (iter->y - clusterY[i]);
          ZSpread += (iter->z - clusterZ[i]) * (iter->z - clusterZ[i]);
          TSpread += (iter->t - clusterTime[i]) * (iter->t - clusterTime[i]);
          TSpreadAll += (iter->t - clusterTimeWeighted[i]) * (iter->t - clusterTimeWeighted[i]);
          float radius = sqrt(pow(iter->x, 2) + pow(iter->y, 2));
          RSpread += pow(radius-sqrt(clusterX[i]*clusterX[i]+clusterY[i]*clusterY[i]),2);


      }
    }

    float a = (m11+m22)/2;
    float b = 0.5*sqrt((m11+m22)*(m11+m22)-4*(m11*m22-m12*m12));
    clusterXSpread.push_back(sqrt(XSpread/(float)clusterSize[i]));
    clusterYSpread.push_back(sqrt(YSpread/(float)clusterSize[i]));
    clusterZSpread.push_back(sqrt(ZSpread/(float)clusterSize[i]));
    clusterRSpread.push_back(sqrt(RSpread/(float)clusterSize[i]));
    clusterDeltaRSpread.push_back(sqrt(DeltaRSpread/(float)clusterSize[i]));
    clusterXYSpread.push_back(sqrt(abs(XYSpread)/(float)clusterSize[i]));
    clusterTimeSpread.push_back(sqrt(TSpread/(float)clusterSize[i]));
    clusterTimeSpreadWeightedAll.push_back(sqrt(TSpreadAll/(float)clusterSize[i]));
    clusterEtaSpread.push_back(sqrt(abs(m11)/clusterSize[i]));
    clusterEtaPhiSpread.push_back(sqrt(abs(m12)/clusterSize[i]));
    clusterPhiSpread.push_back(sqrt(abs(m22)/clusterSize[i]));
    clusterMajorAxis.push_back(sqrt((a+b)/clusterSize[i]));
    clusterMinorAxis.push_back(sqrt((a-b)/clusterSize[i]));

  }

  return 0;
}
//input x,y,z of segments in cluster, avgX, avgY, avgZ
void DBSCAN::sort_clusters()
{
  for(int i = 0; i < nClusters; i++){
    vector<int>dtStations;
    vector<int>dtStations_copy;
    vector<int>dtChambers;
    vector<int>dtChambers_copy;
    vector<int>dtLayersPlus11;
    vector<int>dtLayersPlus12;
    vector<int>dtLayersPlus13;
    vector<int>dtLayersPlus21;
    vector<int>dtLayersPlus22;
    vector<int>dtLayersPlus31;
    vector<int>dtLayersPlus32;
    vector<int>dtLayersPlus41;
    vector<int>dtLayersPlus42;

    vector<int>dtLayersMinus11;
    vector<int>dtLayersMinus12;
    vector<int>dtLayersMinus13;
    vector<int>dtLayersMinus21;
    vector<int>dtLayersMinus22;
    vector<int>dtLayersMinus31;
    vector<int>dtLayersMinus32;
    vector<int>dtLayersMinus41;
    vector<int>dtLayersMinus42;

    vector<int>segment_index;
    int nSegments_Me11 = 0;
    int nSegments_Me12 = 0;

    dtCluster tmpCluster;
    tmpCluster.nDtSegmentChamberPlus11 = 0;
    tmpCluster.nDtSegmentChamberPlus12 = 0;
    tmpCluster.nDtSegmentChamberPlus13 = 0;
    tmpCluster.nDtSegmentChamberPlus21 = 0;
    tmpCluster.nDtSegmentChamberPlus22 = 0;
    tmpCluster.nDtSegmentChamberPlus31 = 0;
    tmpCluster.nDtSegmentChamberPlus32 = 0;
    tmpCluster.nDtSegmentChamberPlus41 = 0;
    tmpCluster.nDtSegmentChamberPlus42 = 0;
    tmpCluster.nDtSegmentChamberMinus11 = 0;
    tmpCluster.nDtSegmentChamberMinus12 = 0;
    tmpCluster.nDtSegmentChamberMinus13 = 0;
    tmpCluster.nDtSegmentChamberMinus21 = 0;
    tmpCluster.nDtSegmentChamberMinus22 = 0;
    tmpCluster.nDtSegmentChamberMinus31 = 0;
    tmpCluster.nDtSegmentChamberMinus32 = 0;
    tmpCluster.nDtSegmentChamberMinus41 = 0;
    tmpCluster.nDtSegmentChamberMinus42 = 0;

    tmpCluster.nDtSegmentStation1 = 0;
    tmpCluster.nDtSegmentStation2 = 0;
    tmpCluster.nDtSegmentStation3 = 0;
    tmpCluster.nDtSegmentStation4 = 0;
    for(unsigned int l=0; l < m_pointSize; l++){
      if (m_points[l].clusterID == i+1){
        dtStations.push_back(m_points[l].station);
        dtChambers.push_back(m_points[l].chamber);
        dtStations_copy.push_back(m_points[l].station);
        dtChambers_copy.push_back(m_points[l].chamber);
        segment_index.push_back(l);

        if (m_points[l].chamber == 11) dtLayersPlus11.push_back(m_points[l].layer);
        if (m_points[l].chamber == 12) dtLayersPlus12.push_back(m_points[l].layer);
        if (m_points[l].chamber == 13) dtLayersPlus13.push_back(m_points[l].layer);
        if (m_points[l].chamber == 21) dtLayersPlus21.push_back(m_points[l].layer);
      	if (m_points[l].chamber == 22) dtLayersPlus22.push_back(m_points[l].layer);
        if (m_points[l].chamber == 31) dtLayersPlus31.push_back(m_points[l].layer);
      	if (m_points[l].chamber == 32) dtLayersPlus32.push_back(m_points[l].layer);
        if (m_points[l].chamber == 41) dtLayersPlus41.push_back(m_points[l].layer);
      	if (m_points[l].chamber == 42) dtLayersPlus42.push_back(m_points[l].layer);

        if (m_points[l].chamber == -11) dtLayersMinus11.push_back(m_points[l].layer);
        if (m_points[l].chamber == -12) dtLayersMinus12.push_back(m_points[l].layer);
        if (m_points[l].chamber == -13) dtLayersMinus13.push_back(m_points[l].layer);
        if (m_points[l].chamber == -21) dtLayersMinus21.push_back(m_points[l].layer);
      	if (m_points[l].chamber == -22) dtLayersMinus22.push_back(m_points[l].layer);
        if (m_points[l].chamber == -31) dtLayersMinus31.push_back(m_points[l].layer);
      	if (m_points[l].chamber == -32) dtLayersMinus32.push_back(m_points[l].layer);
        if (m_points[l].chamber == -41) dtLayersMinus41.push_back(m_points[l].layer);
      	if (m_points[l].chamber == -42) dtLayersMinus42.push_back(m_points[l].layer);

        if (abs(m_points[l].chamber) == 11) nSegments_Me11++;
        if (abs(m_points[l].chamber) == 12) nSegments_Me12++;
      	if (m_points[l].chamber == 11) tmpCluster.nDtSegmentChamberPlus11++;
      	if (m_points[l].chamber == 12) tmpCluster.nDtSegmentChamberPlus12++;
      	if (m_points[l].chamber == 13) tmpCluster.nDtSegmentChamberPlus13++;
      	if (m_points[l].chamber == 21) tmpCluster.nDtSegmentChamberPlus21++;
      	if (m_points[l].chamber == 22) tmpCluster.nDtSegmentChamberPlus22++;
      	if (m_points[l].chamber == 31) tmpCluster.nDtSegmentChamberPlus31++;
      	if (m_points[l].chamber == 32) tmpCluster.nDtSegmentChamberPlus32++;
      	if (m_points[l].chamber == 41) tmpCluster.nDtSegmentChamberPlus41++;
      	if (m_points[l].chamber == 42) tmpCluster.nDtSegmentChamberPlus42++;
      	if (m_points[l].chamber == -11) tmpCluster.nDtSegmentChamberMinus11++;
      	if (m_points[l].chamber == -12) tmpCluster.nDtSegmentChamberMinus12++;
      	if (m_points[l].chamber == -13) tmpCluster.nDtSegmentChamberMinus13++;
      	if (m_points[l].chamber == -21) tmpCluster.nDtSegmentChamberMinus21++;
      	if (m_points[l].chamber == -22) tmpCluster.nDtSegmentChamberMinus22++;
      	if (m_points[l].chamber == -31) tmpCluster.nDtSegmentChamberMinus31++;
      	if (m_points[l].chamber == -32) tmpCluster.nDtSegmentChamberMinus32++;
      	if (m_points[l].chamber == -41) tmpCluster.nDtSegmentChamberMinus41++;
      	if (m_points[l].chamber == -42) tmpCluster.nDtSegmentChamberMinus42++;


        if (m_points[l].station == 1) tmpCluster.nDtSegmentStation1++;
        if (m_points[l].station == 2) tmpCluster.nDtSegmentStation2++;
        if (m_points[l].station == 3) tmpCluster.nDtSegmentStation3++;
        if (m_points[l].station == 4) tmpCluster.nDtSegmentStation4++;

      }
    }
    tmpCluster.x = clusterX[i];
    tmpCluster.y = clusterY[i];
    tmpCluster.z = clusterZ[i];
    tmpCluster.eta = clusterEta[i];
    tmpCluster.phi = clusterPhi[i];
    tmpCluster.t = clusterTime[i];
    tmpCluster.tWeighted = clusterTimeWeighted[i];
    tmpCluster.tTotal = clusterTimeTotal[i];

    tmpCluster.MajorAxis = clusterMajorAxis[i];
    tmpCluster.MinorAxis = clusterMinorAxis[i];
    tmpCluster.XSpread = clusterXSpread[i];
    tmpCluster.XYSpread = clusterXYSpread[i];
    tmpCluster.RSpread = clusterRSpread[i];


    tmpCluster.YSpread = clusterYSpread[i];
    tmpCluster.ZSpread = clusterZSpread[i];
    tmpCluster.TSpread = clusterTimeSpread[i];
    tmpCluster.TSpreadWeighted = clusterTimeSpreadWeighted[i];
    tmpCluster.TSpreadWeightedAll = clusterTimeSpreadWeightedAll[i];

    tmpCluster.EtaPhiSpread = clusterEtaPhiSpread[i];
    tmpCluster.EtaSpread = clusterEtaSpread[i];
    tmpCluster.DeltaRSpread = clusterDeltaRSpread[i];
    tmpCluster.PhiSpread = clusterPhiSpread[i];
    tmpCluster.nDtSegments = clusterSize[i];
    tmpCluster.Me11Ratio = 1.0*nSegments_Me11/clusterSize[i];
    tmpCluster.Me12Ratio = 1.0*nSegments_Me12/clusterSize[i];

    // count number of dtLayers
    std::sort(dtLayersPlus11.begin(), dtLayersPlus11.end());
    std::sort(dtLayersPlus12.begin(), dtLayersPlus12.end());
    std::sort(dtLayersPlus13.begin(), dtLayersPlus13.end());
    std::sort(dtLayersPlus21.begin(), dtLayersPlus21.end());
    std::sort(dtLayersPlus22.begin(), dtLayersPlus22.end());
    std::sort(dtLayersPlus31.begin(), dtLayersPlus31.end());
    std::sort(dtLayersPlus32.begin(), dtLayersPlus32.end());
    std::sort(dtLayersPlus41.begin(), dtLayersPlus41.end());
    std::sort(dtLayersPlus42.begin(), dtLayersPlus42.end());

    std::sort(dtLayersMinus11.begin(), dtLayersMinus11.end());
    std::sort(dtLayersMinus12.begin(), dtLayersMinus12.end());
    std::sort(dtLayersMinus13.begin(), dtLayersMinus13.end());
    std::sort(dtLayersMinus21.begin(), dtLayersMinus21.end());
    std::sort(dtLayersMinus22.begin(), dtLayersMinus22.end());
    std::sort(dtLayersMinus31.begin(), dtLayersMinus31.end());
    std::sort(dtLayersMinus32.begin(), dtLayersMinus32.end());
    std::sort(dtLayersMinus41.begin(), dtLayersMinus41.end());
    std::sort(dtLayersMinus42.begin(), dtLayersMinus42.end());


    tmpCluster.nLayersChamberPlus11 = std::unique(dtLayersPlus11.begin(), dtLayersPlus11.end())-dtLayersPlus11.begin();
    tmpCluster.nLayersChamberPlus12 = std::unique(dtLayersPlus12.begin(), dtLayersPlus12.end())-dtLayersPlus12.begin();
    tmpCluster.nLayersChamberPlus13 = std::unique(dtLayersPlus13.begin(), dtLayersPlus13.end())-dtLayersPlus13.begin();
    tmpCluster.nLayersChamberPlus21 = std::unique(dtLayersPlus21.begin(), dtLayersPlus21.end())-dtLayersPlus21.begin();
    tmpCluster.nLayersChamberPlus22 = std::unique(dtLayersPlus22.begin(), dtLayersPlus22.end())-dtLayersPlus22.begin();
    tmpCluster.nLayersChamberPlus31 = std::unique(dtLayersPlus31.begin(), dtLayersPlus31.end())-dtLayersPlus31.begin();
    tmpCluster.nLayersChamberPlus32 = std::unique(dtLayersPlus32.begin(), dtLayersPlus32.end())-dtLayersPlus32.begin();
    tmpCluster.nLayersChamberPlus41 = std::unique(dtLayersPlus41.begin(), dtLayersPlus41.end())-dtLayersPlus41.begin();
    tmpCluster.nLayersChamberPlus42 = std::unique(dtLayersPlus42.begin(), dtLayersPlus42.end())-dtLayersPlus42.begin();
    tmpCluster.nLayersChamberMinus11 = std::unique(dtLayersMinus11.begin(), dtLayersMinus11.end())-dtLayersMinus11.begin();
    tmpCluster.nLayersChamberMinus12 = std::unique(dtLayersMinus12.begin(), dtLayersMinus12.end())-dtLayersMinus12.begin();
    tmpCluster.nLayersChamberMinus13 = std::unique(dtLayersMinus13.begin(), dtLayersMinus13.end())-dtLayersMinus13.begin();
    tmpCluster.nLayersChamberMinus21 = std::unique(dtLayersMinus21.begin(), dtLayersMinus21.end())-dtLayersMinus21.begin();
    tmpCluster.nLayersChamberMinus22 = std::unique(dtLayersMinus22.begin(), dtLayersMinus22.end())-dtLayersMinus22.begin();
    tmpCluster.nLayersChamberMinus31 = std::unique(dtLayersMinus31.begin(), dtLayersMinus31.end())-dtLayersMinus31.begin();
    tmpCluster.nLayersChamberMinus32 = std::unique(dtLayersMinus32.begin(), dtLayersMinus32.end())-dtLayersMinus32.begin();
    tmpCluster.nLayersChamberMinus41 = std::unique(dtLayersMinus41.begin(), dtLayersMinus41.end())-dtLayersMinus41.begin();
    tmpCluster.nLayersChamberMinus42 = std::unique(dtLayersMinus42.begin(), dtLayersMinus42.end())-dtLayersMinus42.begin();
    tmpCluster.segment_id = segment_index;

    // count the number of chambers and max chamber segments
    std::vector<int>::iterator chamber_it;
    std::sort(dtChambers.begin(), dtChambers.end());
    chamber_it = std::unique(dtChambers.begin(), dtChambers.end());
    dtChambers.resize( std::distance(dtChambers.begin(),chamber_it) );
    int max_chamber = 999; // station with the maximum number of dtsegment in this cluster
    int max_chamber_segment = 0; // station with the maximum number of dtsegment in this cluster
    tmpCluster.nChamber = 0;
    for (unsigned int l = 0; l < dtChambers.size(); l++)
    {
      int counter = 0;
      for(unsigned int j = 0; j < dtChambers_copy.size(); j ++)
      {
        if (dtChambers_copy[j] == dtChambers[l]) counter++;
      }
      if (counter>max_chamber_segment)
      {
        max_chamber_segment = counter;
        max_chamber = dtChambers[l];
      }
      if(counter>5)tmpCluster.nChamber++;
    }
    tmpCluster.maxChamber = max_chamber;
    tmpCluster.maxChamberSegment = max_chamber_segment;
    // count the number of chambers and max chamber segments
    std::vector<int>::iterator station_it;
    std::sort(dtStations.begin(), dtStations.end());
    station_it = std::unique(dtStations.begin(), dtStations.end());
    dtStations.resize( std::distance(dtStations.begin(),station_it) );//list of unique stations
    int max_station = 999; // station with the maximum number of dtsegment in this cluster
    int max_station_segment = 0; // station with the maximum number of dtsegment in this cluster
    tmpCluster.nStation = 0;
    tmpCluster.nStation5 = 0;
    tmpCluster.nStation10 = 0;
    tmpCluster.nStation10perc = 0;
    tmpCluster.avgStation = 0.0;
    tmpCluster.avgStation5 = 0.0;
    tmpCluster.avgStation10 = 0.0;
    tmpCluster.avgStation10perc = 0.0;
    int nSeg10perc = 0;
    int nSeg5 = 0;
    int nSeg10 = 0;
    for (unsigned int l = 0; l < dtStations.size(); l++)
    {
      int counter = 0;
      for(unsigned int j = 0; j < dtStations_copy.size(); j ++)
      {
        if (dtStations_copy[j] == dtStations[l]) counter++;
      }
      if (counter>max_station_segment)
      {
        max_station_segment = counter;
        max_station = dtStations[l];
      }
      tmpCluster.nStation++;
      tmpCluster.avgStation += counter * dtStations[l];

      if(counter>=5.0){
        tmpCluster.avgStation5 += counter * dtStations[l];
        tmpCluster.nStation5++;
        nSeg5 += counter;

      }
      if(counter>=10.0){
        tmpCluster.avgStation10 += counter * dtStations[l];
        tmpCluster.nStation10++;
        nSeg10 += counter;
      }
      if(1.0*counter/clusterSize[i] > 0.1)
      {
        tmpCluster.avgStation10perc += counter * dtStations[l];
        tmpCluster.nStation10perc++;
        nSeg10perc += counter;
      }
    }
    tmpCluster.avgStation10perc = 1.0* tmpCluster.avgStation10perc/nSeg10perc;
    tmpCluster.avgStation5 = 1.0* tmpCluster.avgStation5/nSeg5;
    tmpCluster.avgStation10 = 1.0* tmpCluster.avgStation10/nSeg10;
    tmpCluster.avgStation = 1.0* tmpCluster.avgStation/clusterSize[i];

    tmpCluster.maxStation = max_station;
    tmpCluster.maxStationSegment = max_station_segment;

    DtCluster.push_back(tmpCluster);


  }

  //sort the clusters by size
  sort(DtCluster.begin(), DtCluster.end(), largest_nDt_cluster);

}


void DBSCAN::merge_clusters()
{
  // clear all the cluster variables
  //change cluster ID of points
  bool modified = true;
  while(modified){
    modified = false;
    float mindR = 15;
    int cluster1 = 999;
    int cluster2 = 999;

    for(unsigned int i = 0; i < clusterEta.size(); i++){
      for(unsigned int j = i+1; j < clusterEta.size(); j++){
        float current_dR = deltaR(clusterEta[i], clusterPhi[i], clusterEta[j], clusterPhi[j]);
        if(current_dR<mindR)
        {
          mindR = current_dR;
          cluster1 = i;
          cluster2 = j;

        }
      }
    }
    if (mindR < MERGE_CLUSTER_DR){
      vector<Point>::iterator iter;
      int count = 0;
      for(iter = m_points.begin(); iter != m_points.end(); ++iter)
      {
        if ( iter->clusterID == cluster2+1 ){
          iter->clusterID = cluster1+1;
          count++;
        }
        if ( iter->clusterID > cluster2+1 )iter->clusterID = iter->clusterID-1;
      }
      clusterEta.erase(clusterEta.begin() + cluster2);
      clusterPhi.erase(clusterPhi.begin() + cluster2);
      nClusters--;
      modified = true;
      // can't use dtCluster, because its sorted, but the other vectors and clusterID are not sorted.
    }
  }
  clear_clusters(); // clear all the vectors, but nClusters is kept to keep track of the number of clusters.
}

double DBSCAN::deltaPhi(double phi1, double phi2)
{
  double dphi = phi1-phi2;
  while (dphi > TMath::Pi())
  {
    dphi -= TMath::TwoPi();
  }
  while (dphi <= -TMath::Pi())
  {
    dphi += TMath::TwoPi();
  }
  return dphi;
};
double DBSCAN::deltaR(double eta1, double phi1, double eta2, double phi2) {
  double dphi = deltaPhi(phi1,phi2);
  double deta = eta1 - eta2;
  return sqrt( dphi*dphi + deta*deta);
}
