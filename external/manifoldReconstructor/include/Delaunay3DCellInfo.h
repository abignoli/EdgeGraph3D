/*
 * [Imported from external project]
 *
 * OpenMvgParser.cpp
 *
 *		Source: Manifold Reconstructor [https://github.com/andresax/Manifold-Reconstructor]
 *      Author: Andrea Romanoni
 */

#ifndef DELAUNAY3DCELLINFO_H_
#define DELAUNAY3DCELLINFO_H_

#define NO_HEURISTIC_K
#define HEURISTIC_K 5
//#define HEURISTIC_K 1


#include <Eigen/Core>
#include <set>
#include <vector>
#include <glm.hpp>



class Delaunay3DCellInfo {
public:
  // Sorted-Set related struct
  struct FSConstraint;
  /*struct LtConstraint {
    bool operator()(const std::pair<int, int> x, const std::pair<int, int> y) const {
      if (x.first < y.first)
        return true;
      else
        return (x.first == y.first) && (x.second < y.second);
    }
  };*/
  struct LtFSConstraint {
    bool operator()(const FSConstraint x, const FSConstraint y) const {
      if (x.first < y.first)
        return true;
      else
        return (x.first == y.first) && (x.vertexIdx < y.vertexIdx);
    }
  };
  struct FSConstraint {
    FSConstraint() {
      first = -1;
      second = -1;
      vote = -1;
      vertexIdx= second;
      fNearestNeighborDist = std::numeric_limits<float>::infinity();
    }
    FSConstraint(int camIndex, int featureIndex) {
      first = camIndex;
      second = featureIndex;
      vote = 1.0;
      vertexIdx= second;
      fNearestNeighborDist = std::numeric_limits<float>::infinity();
    }
    FSConstraint(int camIndex, int featureIndex, float voteVal) {
      first = camIndex;
      second = featureIndex;
      vertexIdx= second;
      vote = voteVal;
      fNearestNeighborDist = std::numeric_limits<float>::infinity();
    }
    FSConstraint(int camIndex, int featureIndex, int vertidx,float voteVal) {
      first = camIndex;
      second = featureIndex;
      vertexIdx = vertidx;
      vote = voteVal;
      fNearestNeighborDist = std::numeric_limits<float>::infinity();
    }
    FSConstraint(const FSConstraint & ref) {
      first = ref.first;
      second = ref.second;
      vote = ref.vote;
      vertexIdx = ref.vertexIdx;
      pNearestNeighbor = ref.pNearestNeighbor;
      fNearestNeighborDist = ref.fNearestNeighborDist;
    }
    FSConstraint & operator=(const FSConstraint & ref) {
      if (&ref != this) {
        first = ref.first;
        second = ref.second;
        vote = ref.vote;
        pNearestNeighbor = ref.pNearestNeighbor;
        fNearestNeighborDist = ref.fNearestNeighborDist;
      }
      return *this;
    }
    operator std::pair<int, int>() const {
      return std::make_pair(first, second);
    }

    int first;
    int second;
    int vertexIdx;
    float vote;
    mutable std::set<FSConstraint, LtFSConstraint>::iterator pNearestNeighbor;
    mutable float fNearestNeighborDist;

    // Trick c++ "const"-ness so that we can change the nearest neighbor information in the set from an iterator:
    void setNearestNeighbor(const std::set<FSConstraint, LtFSConstraint>::iterator itNearest, const float dist) const {
      pNearestNeighbor = itNearest;
      fNearestNeighborDist = dist;
    }
    void resetNearestNeighborDist() const {
      fNearestNeighborDist = std::numeric_limits<float>::infinity();
    }
  };

  // Constructors (it must be default-constructable)
  Delaunay3DCellInfo();
  Delaunay3DCellInfo(const Delaunay3DCellInfo & ref);


  virtual ~Delaunay3DCellInfo();

  // Getters
  int getVoteCount() const {
    return m_voteCount;
  }
  float getVoteCountProb() const {
    return m_voteCountProb;
  }

  bool isBoundary() const {
    return boundary;
  }
  bool iskeptManifold() const {
    return keepManifold;
  }

  bool isTemporaryInside() const {
    return temporary_Inside_;
  }
  bool isToBeTested() const {
    return toBeTested_;
  }
  bool isShrinked() const {
    return shrinked_;
  }
  bool isGalpha() const {
    return Galpha_;
  }
  const std::set<FSConstraint, LtFSConstraint> & getIntersections() const {
    return m_setIntersections;
  }
  bool isNew() const {
    return m_bNew;
  }

  // Setters
  void setVoteCount(const int voteCount) {
    m_voteCount = voteCount;
  }

  void setBoundary(bool value) {
    boundary = value;
  }
  void setKeptManifold(bool value) {
    keepManifold = value;
  }
  void setShrinked(bool value) {
    shrinked_ = value;
  }
  void setGalpha(bool value) {
    Galpha_ = value;
  }
  void setTemporaryInside(bool value) {
    temporary_Inside_ = value;
  }
  void setToBeTested(bool value) {
    toBeTested_ = value;
  }
  void setVoteCountProb(const float voteCountProb) {
    m_voteCountProb = voteCountProb;
  }
  void setIntersections(const std::set<FSConstraint, LtFSConstraint> & ref) {
    m_setIntersections = ref;
  }

  void setIdxManifold(bool value, int i) {
    idxFacetNotManifold[i] = value;
  }

  bool getIdxManifold(int i) {
    return idxFacetNotManifold[i];
  }

  // Public Methods
  void incrementVoteCount() {
    m_voteCount++;
  }

  void incrementVoteCount(int num) {
    m_voteCount += num;
  }

  void incrementVoteCountProb(float incr) {
    m_voteCountProb += incr;
  }
  void decrementVoteCount() {
    if (m_voteCount > 0)
      m_voteCount--;
  }
  void decrementVoteCount(int num) {
    if (m_voteCount - num > 0)
      m_voteCount -= num;
    else
      m_voteCount = 0;
  }
  void decrementVoteCountProb(float incr) {
    if (m_voteCountProb > incr)
      m_voteCountProb -= incr;
  }
  bool isKeptByVoteCount(const int nVoteThresh = 1) const {
    if (getVoteCount() < nVoteThresh)
      return true;
    return false;
  }
  bool isKeptByVoteCountProb(const float nVoteThresh = 1.0) const {
    if (getVoteCountProb() < nVoteThresh)
      return true;
    return false;
  }

  /* template<class T>
   void addIntersection(int camIndex, int featureIndex, const vector<T> & vecVertexHandles, const vector<Matrix> & vecCamCenters) {
   m_setIntersections.insert(m_setIntersections.end(), FSConstraint(camIndex, featureIndex));
   }
   template<class T>
   void addIntersection(int camIndex, int featureIndex, float vote, const vector<T> & vecVertexHandles, const vector<Matrix> & vecCamCenters) {
   m_setIntersections.insert(m_setIntersections.end(), FSConstraint(camIndex, featureIndex,vote));
   }*/
#ifdef NO_HEURISTIC_K
  template<class T>
  void addIntersection(int camIndex, int featureIndex, float vote, const std::vector<T> & vecVertexHandles, const std::vector<glm::vec3> & vecCamCenters) {
    m_setIntersections.insert(m_setIntersections.end(), FSConstraint(camIndex, featureIndex, vote));
  }


  template<class T>
  void addIntersection(int camIndex, int featureIndex, float vote, int vertexidx, const std::vector<T> & vecVertexHandles, const std::vector<glm::vec3> & vecCamCenters) {
    m_setIntersections.insert(m_setIntersections.end(), FSConstraint(camIndex, featureIndex, vertexidx, vote));
  }

  template<class T>
  void addIntersection(int camIndex, int featureIndex, const std::vector<T> & vecVertexHandles, const std::vector<glm::vec3> & vecCamCenters) {
    m_setIntersections.insert(m_setIntersections.end(), FSConstraint(camIndex, featureIndex));
  }
#else
  template<class T>
  void addIntersection(int camIndex, int featureIndex, const std::vector<T> & vecVertexHandles, const std::vector<glm::vec3> & vecCamCenters) {
    FSConstraint incoming(camIndex, featureIndex);
    if ((int) m_setIntersections.size() < HEURISTIC_K) {
      // The constraint set is not full, so insert the incoming free-space constraint and update nearest neighbor info.
      std::set<FSConstraint, LtFSConstraint>::iterator it, itIncoming;
      itIncoming = m_setIntersections.insert(m_setIntersections.end(), incoming);
      for (it = m_setIntersections.begin(); it != m_setIntersections.end(); it++) {
        if (it == itIncoming)
          continue;
        // Asymmetric metric:
        float curDist = distFSConstraint(*itIncoming, *it, vecVertexHandles, vecCamCenters);
        float curDist2 = distFSConstraint(*it, *itIncoming, vecVertexHandles, vecCamCenters);
        // Update incoming
        if (curDist < itIncoming->fNearestNeighborDist)
          itIncoming->setNearestNeighbor(it, curDist);
        // Update *it:
        if (curDist2 < it->fNearestNeighborDist)
          it->setNearestNeighbor(itIncoming, curDist2);
      }
    } else {
      // The constraint set is full, so apply the spatial cover heuristic to determine whether or not to insert the incoming free-space constraint

      // Quick return / rejection on the case that m_nMaxConstraintsKept == 1.
      if (HEURISTIC_K == 1)
        return;

      float minDist = std::numeric_limits<float>::infinity();
      std::set<FSConstraint, LtFSConstraint>::iterator it, it2, itEject;
      for (it = m_setIntersections.begin(); it != m_setIntersections.end(); it++) {
        float curDist = distFSConstraint(incoming, *it, vecVertexHandles, vecCamCenters);
        if (curDist < it->fNearestNeighborDist)
          break; // REJECT
        float curDist2 = distFSConstraint(*it, incoming, vecVertexHandles, vecCamCenters);
        if (curDist2 < it->fNearestNeighborDist)
          break; // REJECT
        // Update incoming
        if (curDist < incoming.fNearestNeighborDist)
          incoming.setNearestNeighbor(it, curDist);
        // Update minDist & itEject
        if (it->fNearestNeighborDist < minDist) {
          minDist = it->fNearestNeighborDist;
          itEject = it;
        }
      }

      if (it == m_setIntersections.end()) {
        // No rejection, so insert incoming and evict itEject.

        // For an asymmetric metric, incoming might have its nearest neighbor ejected.  If so, compute the 2nd nearest neighbor
        if (incoming.pNearestNeighbor == itEject) {
          incoming.fNearestNeighborDist = std::numeric_limits<float>::infinity();
          for (it2 = m_setIntersections.begin(); it2 != m_setIntersections.end(); it2++) {
            if (it2 == itEject)
              continue;
            float curDist = distFSConstraint(incoming, *it2, vecVertexHandles, vecCamCenters);
            if (curDist < incoming.fNearestNeighborDist)
              incoming.setNearestNeighbor(it2, curDist);
          }
        }

        // Recompute nearest neighbors that previously pointed to itEject.
        for (it2 = m_setIntersections.begin(); it2 != m_setIntersections.end(); it2++) {
          if (it2->pNearestNeighbor == itEject) { // implicity "continue;"'s if it2 == itEject
            // Recompute the nearest neighbor for it2:
            it2->resetNearestNeighborDist();
            for (it = m_setIntersections.begin(); it != m_setIntersections.end(); it++) {
              if (it == itEject || it == it2)
                continue;
              float curDist = distFSConstraint(*it2, *it, vecVertexHandles, vecCamCenters);
              if (curDist < it2->fNearestNeighborDist)
                it2->setNearestNeighbor(it, curDist);
            }
          }
        }

        // Finally erase itEject and insert incoming
        m_setIntersections.erase(itEject);
        m_setIntersections.insert(m_setIntersections.end(), incoming);
      }
    }
  }

  template<class T>
  void addIntersection(int camIndex, int featureIndex, float vote, const std::vector<T> & vecVertexHandles, const std::vector<glm::vec3> & vecCamCenters) {
    FSConstraint incoming(camIndex, featureIndex, vote);
    if ((int) m_setIntersections.size() < HEURISTIC_K) {
      // The constraint set is not full, so insert the incoming free-space constraint and update nearest neighbor info.
      std::set<FSConstraint, LtFSConstraint>::iterator it, itIncoming;
      itIncoming = m_setIntersections.insert(m_setIntersections.end(), incoming);
      for (it = m_setIntersections.begin(); it != m_setIntersections.end(); it++) {
        if (it == itIncoming)
          continue;
        // Asymmetric metric:
        float curDist = distFSConstraint(*itIncoming, *it, vecVertexHandles, vecCamCenters);
        float curDist2 = distFSConstraint(*it, *itIncoming, vecVertexHandles, vecCamCenters);
        // Update incoming
        if (curDist < itIncoming->fNearestNeighborDist)
          itIncoming->setNearestNeighbor(it, curDist);
        // Update *it:
        if (curDist2 < it->fNearestNeighborDist)
          it->setNearestNeighbor(itIncoming, curDist2);
      }
    } else {
      // The constraint set is full, so apply the spatial cover heuristic to determine whether or not to insert the incoming free-space constraint

      // Quick return / rejection on the case that m_nMaxConstraintsKept == 1.
      if (HEURISTIC_K == 1)
        return;

      float minDist = std::numeric_limits<float>::infinity();
      std::set<FSConstraint, LtFSConstraint>::iterator it, it2, itEject;
      for (it = m_setIntersections.begin(); it != m_setIntersections.end(); it++) {
        float curDist = distFSConstraint(incoming, *it, vecVertexHandles, vecCamCenters);
        if (curDist < it->fNearestNeighborDist)
          break; // REJECT
        float curDist2 = distFSConstraint(*it, incoming, vecVertexHandles, vecCamCenters);
        if (curDist2 < it->fNearestNeighborDist)
          break; // REJECT
        // Update incoming
        if (curDist < incoming.fNearestNeighborDist)
          incoming.setNearestNeighbor(it, curDist);
        // Update minDist & itEject
        if (it->fNearestNeighborDist < minDist) {
          minDist = it->fNearestNeighborDist;
          itEject = it;
        }
      }

      if (it == m_setIntersections.end()) {
        // No rejection, so insert incoming and evict itEject.

        // For an asymmetric metric, incoming might have its nearest neighbor ejected.  If so, compute the 2nd nearest neighbor
        if (incoming.pNearestNeighbor == itEject) {
          incoming.fNearestNeighborDist = std::numeric_limits<float>::infinity();
          for (it2 = m_setIntersections.begin(); it2 != m_setIntersections.end(); it2++) {
            if (it2 == itEject)
              continue;
            float curDist = distFSConstraint(incoming, *it2, vecVertexHandles, vecCamCenters);
            if (curDist < incoming.fNearestNeighborDist)
              incoming.setNearestNeighbor(it2, curDist);
          }
        }

        // Recompute nearest neighbors that previously pointed to itEject.
        for (it2 = m_setIntersections.begin(); it2 != m_setIntersections.end(); it2++) {
          if (it2->pNearestNeighbor == itEject) { // implicity "continue;"'s if it2 == itEject
            // Recompute the nearest neighbor for it2:
            it2->resetNearestNeighborDist();
            for (it = m_setIntersections.begin(); it != m_setIntersections.end(); it++) {
              if (it == itEject || it == it2)
                continue;
              float curDist = distFSConstraint(*it2, *it, vecVertexHandles, vecCamCenters);
              if (curDist < it2->fNearestNeighborDist)
                it2->setNearestNeighbor(it, curDist);
            }
          }
        }

        // Finally erase itEject and insert incoming
        m_setIntersections.erase(itEject);
        m_setIntersections.insert(m_setIntersections.end(), incoming);
      }
    }
  }
#endif
#ifdef NO_HEURISTIC_K
  template<class T>
  void removeIntersection(int camIndex, int featureIndex, const std::vector<T> & vecVertexHandles, const std::vector<glm::vec3> & vecCamCenters) {
    m_setIntersections.erase(FSConstraint(camIndex, featureIndex));
  }

  template<class T>
  void removeIntersection(int camIndex, int featureIndex, float vote, const std::vector<T> & vecVertexHandles,
      const std::vector<glm::vec3> & vecCamCenters) {
    m_setIntersections.erase(FSConstraint(camIndex, featureIndex, vote));
  }
#else
  template<class T>
  void removeIntersection(int camIndex, int featureIndex, const std::vector<T> & vecVertexHandles, const std::vector<glm::vec3> & vecCamCenters) {
    if ((int) m_setIntersections.size() <= 1) {
      // No nearest neighbor info needs to be updated
      m_setIntersections.erase(FSConstraint(camIndex, featureIndex));
    } else {
      // The nearest neighbor info needs to be updated
      std::set<FSConstraint, LtFSConstraint>::iterator it, it2, itEject;

      itEject = m_setIntersections.find(FSConstraint(camIndex, featureIndex));
      if (itEject == m_setIntersections.end())
        return;                   // wasn't in the set to begin with

      for (it = m_setIntersections.begin(); it != m_setIntersections.end(); it++) {
        if (it == itEject)
          continue;
        if (it->pNearestNeighbor == itEject) {
          // Then recompute the nearest neighbor for it:
          it->resetNearestNeighborDist();
          for (it2 = m_setIntersections.begin(); it2 != m_setIntersections.end(); it2++) {
            if (it2 == itEject || it2 == it)
              continue;
            float curDist = distFSConstraint(*it, *it2, vecVertexHandles, vecCamCenters);
            if (curDist < it->fNearestNeighborDist)
              it->setNearestNeighbor(it2, curDist);
          }
        }
      }

      // Finally, erase it.
      m_setIntersections.erase(itEject);
    }
  }
  template<class T>
  void removeIntersection(int camIndex, int featureIndex, float vote, const std::vector<T> & vecVertexHandles, const std::vector<glm::vec3> & vecCamCenters) {
    if ((int) m_setIntersections.size() <= 1) {
      // No nearest neighbor info needs to be updated
      m_setIntersections.erase(FSConstraint(camIndex, featureIndex, vote));
    } else {
      // The nearest neighbor info needs to be updated
      std::set<FSConstraint, LtFSConstraint>::iterator it, it2, itEject;

      itEject = m_setIntersections.find(FSConstraint(camIndex, featureIndex, vote));
      if (itEject == m_setIntersections.end())
        return;                   // wasn't in the set to begin with

      for (it = m_setIntersections.begin(); it != m_setIntersections.end(); it++) {
        if (it == itEject)
          continue;
        if (it->pNearestNeighbor == itEject) {
          // Then recompute the nearest neighbor for it:
          it->resetNearestNeighborDist();
          for (it2 = m_setIntersections.begin(); it2 != m_setIntersections.end(); it2++) {
            if (it2 == itEject || it2 == it)
              continue;
            float curDist = distFSConstraint(*it, *it2, vecVertexHandles, vecCamCenters);
            if (curDist < it->fNearestNeighborDist)
              it->setNearestNeighbor(it2, curDist);
          }
        }
      }

      // Finally, erase it.
      m_setIntersections.erase(itEject);
    }
  }
#endif

  int getNumIntersection() {
    return m_setIntersections.size();
  }

  /*  template<class T>
   void removeIntersection(int camIndex, int featureIndex, const vector<T> & vecVertexHandles, const vector<Matrix> & vecCamCenters) {
   m_setIntersections.erase(FSConstraint(camIndex, featureIndex));
   }

   template<class T>
   void removeIntersection(int camIndex, int featureIndex, float vote, const vector<T> & vecVertexHandles,
   const vector<Matrix> & vecCamCenters) {
   m_setIntersections.erase(FSConstraint(camIndex, featureIndex, vote));
   }*/

  template<class T>
  float distFSConstraint(const FSConstraint & x, const FSConstraint & y, const std::vector<T> & vecVertexHandles, const std::vector<glm::vec3> & vecCamCenters) {
    return distFSConstraintTriangleAreaAaron(x, y, vecVertexHandles, vecCamCenters);
  }
  void clearIntersections() {
    m_setIntersections.clear();
  }
  void markOld() {
    m_bNew = false;
  }
  void markNew() {
    m_bNew = true;
  }

  // Operators (It must be assignable)
  Delaunay3DCellInfo & operator=(const Delaunay3DCellInfo & rhs) {
    if (this != &rhs) {
      setVoteCount(rhs.getVoteCount());
      setIntersections(rhs.getIntersections());
      if (!rhs.isNew())
        markOld();
    }
    return *this;
  }

private:
  // Private Methods
  template<class T>
  float distFSConstraintTriangleAreaAaron(const FSConstraint & x, const FSConstraint & y, const std::vector<T> & vecVertexHandles,
      const std::vector<glm::vec3> & vecCamCenters) {
    // Asymmetric distance heuristic.
    // Sum of two triangle areas, use the base segment PQ as constraint x, and the two points from y as R1 and R2.
    // Note: For efficiency, to avoid unnecessary division by 2 and square-roots, use the sum of twice-the-areas squared = squared area of parallelograms.
    const Eigen::Vector3f & P = Eigen::Vector3f(vecCamCenters[x.first].x,vecCamCenters[x.first].y,
        vecCamCenters[x.first].z);
    Eigen::Vector3f Q;
    Q(0) = vecVertexHandles[x.second].position.x();
    Q(1) = vecVertexHandles[x.second].position.y();
    Q(2) = vecVertexHandles[x.second].position.z();
    const Eigen::Vector3f & R1 = Eigen::Vector3f(vecCamCenters[y.first].x,vecCamCenters[y.first].y,
        vecCamCenters[y.first].z);
    Eigen::Vector3f R2(3, 1);
    R2(0) = vecVertexHandles[y.second].position.x();
    R2(1) = vecVertexHandles[y.second].position.y();
    R2(2) = vecVertexHandles[y.second].position.z();

    // Vector distances
    Eigen::Vector3f PQ(Q - P);
    Eigen::Vector3f PR1(R1 - P);
    Eigen::Vector3f PR2(R2 - P);

    // Sum of squared areas of parallelograms
    Eigen::Vector3f PQxPR1(PQ.cross(PR1));
    Eigen::Vector3f PQxPR2(PQ.cross(PR2));
    return PQxPR1.dot(PQxPR1) + PQxPR2.dot(PQxPR2);
  }

  // Private Members
  int m_voteCount;
  bool boundary;
  bool keepManifold;
  bool toBeTested_;
  float m_voteCountProb;
  std::set<FSConstraint, LtFSConstraint> m_setIntersections;

  bool m_bNew;
  bool idxFacetNotManifold[4];
  bool shrinked_;
  bool Galpha_;
  bool temporary_Inside_;
};

#endif /* DELAUNAY3DCELLINFO_H_ */
