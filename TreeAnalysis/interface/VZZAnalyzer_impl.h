#ifndef VZZANALYZER_IMPL_H
#define VZZANALYZER_IMPL_H

/** \Implementation of template functions of VZZAnalyzer
 *
 *  $Date: 2020/08/27 12:00:00 $
 *  $Revision: 1.0 $
 *  \author A. Mecca alberto.mecca@edu.unito.it
 */

#include "VVXAnalysis/TreeAnalysis/interface/VZZAnalyzer.h"
#include "VVXAnalysis/Commons/interface/Utilities.h"

#include <utility>  // std::pair, std::make_pair()

using std::vector;
using std::pair;
using std::make_pair;
using namespace phys;


template <class P, class R = Boson<Particle>> // P = Jet or Particle
const P* VZZAnalyzer::closestSing(vector<P>* cands, const R& reference, size_t& k){
	if(cands->size() == 0) return nullptr;
	if(cands->size() == 1) { k = 0; return &(cands->front());}
	
	//std::sort( cands->begin(), cands->end(), phys::DeltaRComparator(reference) );
	float minDR = 0.4; //starting value = the threshold we use
	for(size_t i = 0; i < cands->size(); ++i){
	float DR = physmath::deltaR( cands->at(i), reference );
		if(DR < minDR){
			minDR = DR;
			k = i;
		}
	}
	if(k < 99 && physmath::deltaR(reference, cands->at(k)) < 0.4) //cands->front()) < 0.4 )
		return &(cands->at(k));
	else return nullptr;
}


template <class P, class R = Boson<Particle>> // P = Jet or Particle
Boson<P>* VZZAnalyzer::closestPair(vector<P>* cands, const R& reference){
	if(cands->size() < 2) return nullptr;
	else{
		if(cands->size() == 2){
			Boson<P>* res = new Boson<P>( cands->at(0), cands->at(1) );
			if(physmath::deltaR( *res, reference ) < 0.4 )
				return res;
			else return nullptr;
		}
		else{
			//Find the pair with the closest p4
			pair<size_t, size_t> indices(0,0);
			float minDR = 0.4; //starting value = the threshold we use
			for(size_t i = 0; i < cands->size(); ++i)
				for(size_t j = i+1; j < cands->size(); ++j){
					TLorentzVector p4Cand = cands->at(i).p4() + cands->at(j).p4();
					float DR = physmath::deltaR( p4Cand, reference.p4() );
					if(DR < minDR){
						minDR = DR;
						indices = std::make_pair(i,j);
					}
				}
			if(indices.second != 0) //then we've found a good pair
				return new Boson<P>( cands->at(indices.first), cands->at(indices.second) );
			else return nullptr;
		}
	}
}

#endif
