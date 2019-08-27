// Do NOT change. Changes will be lost next time file is generated

#define R__DICTIONARY_FILENAME G__TreeAnalysis

/*******************************************************************/
#include <stddef.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#define G__DICTIONARY
#include "RConfig.h"
#include "TClass.h"
#include "TDictAttributeMap.h"
#include "TInterpreter.h"
#include "TROOT.h"
#include "TBuffer.h"
#include "TMemberInspector.h"
#include "TInterpreter.h"
#include "TVirtualMutex.h"
#include "TError.h"

#ifndef G__ROOT
#define G__ROOT
#endif

#include "RtypesImp.h"
#include "TIsAProxy.h"
#include "TFileMergeInfo.h"
#include <algorithm>
#include "TCollectionProxyInfo.h"
/*******************************************************************/

#include "TDataMember.h"

// Since CINT ignores the std namespace, we need to do so in this file.
namespace std {} using namespace std;

// Header files passed as explicit arguments
#include "/home/cristiano/VVXAnalysis/TreeAnalysis/../..//VVXAnalysis/DataFormats/interface/DiBoson.h"
#include "/home/cristiano/VVXAnalysis/TreeAnalysis/../..//VVXAnalysis/DataFormats/interface/TypeDefs.h"
#include "/home/cristiano/VVXAnalysis/TreeAnalysis/../..//VVXAnalysis/DataFormats/interface/GenEventWeights.h"
#include "/home/cristiano/VVXAnalysis/TreeAnalysis/../..//VVXAnalysis/DataFormats/interface/MELA.h"
#include "/home/cristiano/VVXAnalysis/TreeAnalysis/../..//VVXAnalysis/DataFormats/interface/Jet.h"
#include "/home/cristiano/VVXAnalysis/TreeAnalysis/../..//VVXAnalysis/DataFormats/interface/Boson.h"
#include "/home/cristiano/VVXAnalysis/TreeAnalysis/../..//VVXAnalysis/DataFormats/interface/GenStatusBit.h"
#include "/home/cristiano/VVXAnalysis/TreeAnalysis/../..//VVXAnalysis/DataFormats/interface/Lepton.h"
#include "/home/cristiano/VVXAnalysis/TreeAnalysis/../..//VVXAnalysis/DataFormats/interface/Electron.h"
#include "/home/cristiano/VVXAnalysis/TreeAnalysis/../..//VVXAnalysis/DataFormats/interface/Particle.h"
#include "VVXAnalysis/DataFormats/interface/Particle.h"
#include "VVXAnalysis/DataFormats/interface/Lepton.h"
#include "VVXAnalysis/DataFormats/interface/Jet.h"
#include "VVXAnalysis/DataFormats/interface/Electron.h"
#include "VVXAnalysis/DataFormats/interface/Boson.h"
#include "VVXAnalysis/DataFormats/interface/DiBoson.h"
#include "VVXAnalysis/DataFormats/interface/GenEventWeights.h"
#include "VVXAnalysis/DataFormats/interface/MELA.h"

// Header files passed via #pragma extra_include

namespace ROOT {
   static TClass *pairlEphyscLcLBosonlEphyscLcLLeptongRcOphyscLcLLeptongR_Dictionary();
   static void pairlEphyscLcLBosonlEphyscLcLLeptongRcOphyscLcLLeptongR_TClassManip(TClass*);
   static void *new_pairlEphyscLcLBosonlEphyscLcLLeptongRcOphyscLcLLeptongR(void *p = 0);
   static void *newArray_pairlEphyscLcLBosonlEphyscLcLLeptongRcOphyscLcLLeptongR(Long_t size, void *p);
   static void delete_pairlEphyscLcLBosonlEphyscLcLLeptongRcOphyscLcLLeptongR(void *p);
   static void deleteArray_pairlEphyscLcLBosonlEphyscLcLLeptongRcOphyscLcLLeptongR(void *p);
   static void destruct_pairlEphyscLcLBosonlEphyscLcLLeptongRcOphyscLcLLeptongR(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const pair<phys::Boson<phys::Lepton>,phys::Lepton>*)
   {
      pair<phys::Boson<phys::Lepton>,phys::Lepton> *ptr = 0;
      static ::TVirtualIsAProxy* isa_proxy = new ::TIsAProxy(typeid(pair<phys::Boson<phys::Lepton>,phys::Lepton>));
      static ::ROOT::TGenericClassInfo 
         instance("pair<phys::Boson<phys::Lepton>,phys::Lepton>", "string", 96,
                  typeid(pair<phys::Boson<phys::Lepton>,phys::Lepton>), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &pairlEphyscLcLBosonlEphyscLcLLeptongRcOphyscLcLLeptongR_Dictionary, isa_proxy, 4,
                  sizeof(pair<phys::Boson<phys::Lepton>,phys::Lepton>) );
      instance.SetNew(&new_pairlEphyscLcLBosonlEphyscLcLLeptongRcOphyscLcLLeptongR);
      instance.SetNewArray(&newArray_pairlEphyscLcLBosonlEphyscLcLLeptongRcOphyscLcLLeptongR);
      instance.SetDelete(&delete_pairlEphyscLcLBosonlEphyscLcLLeptongRcOphyscLcLLeptongR);
      instance.SetDeleteArray(&deleteArray_pairlEphyscLcLBosonlEphyscLcLLeptongRcOphyscLcLLeptongR);
      instance.SetDestructor(&destruct_pairlEphyscLcLBosonlEphyscLcLLeptongRcOphyscLcLLeptongR);
      return &instance;
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal((const pair<phys::Boson<phys::Lepton>,phys::Lepton>*)0x0); R__UseDummy(_R__UNIQUE_DICT_(Init));

   // Dictionary for non-ClassDef classes
   static TClass *pairlEphyscLcLBosonlEphyscLcLLeptongRcOphyscLcLLeptongR_Dictionary() {
      TClass* theClass =::ROOT::GenerateInitInstanceLocal((const pair<phys::Boson<phys::Lepton>,phys::Lepton>*)0x0)->GetClass();
      pairlEphyscLcLBosonlEphyscLcLLeptongRcOphyscLcLLeptongR_TClassManip(theClass);
   return theClass;
   }

   static void pairlEphyscLcLBosonlEphyscLcLLeptongRcOphyscLcLLeptongR_TClassManip(TClass* ){
   }

} // end of namespace ROOT

namespace ROOT {
   static void *new_physcLcLParticle(void *p = 0);
   static void *newArray_physcLcLParticle(Long_t size, void *p);
   static void delete_physcLcLParticle(void *p);
   static void deleteArray_physcLcLParticle(void *p);
   static void destruct_physcLcLParticle(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const ::phys::Particle*)
   {
      ::phys::Particle *ptr = 0;
      static ::TVirtualIsAProxy* isa_proxy = new ::TInstrumentedIsAProxy< ::phys::Particle >(0);
      static ::ROOT::TGenericClassInfo 
         instance("phys::Particle", ::phys::Particle::Class_Version(), "Particle.h", 29,
                  typeid(::phys::Particle), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &::phys::Particle::Dictionary, isa_proxy, 4,
                  sizeof(::phys::Particle) );
      instance.SetNew(&new_physcLcLParticle);
      instance.SetNewArray(&newArray_physcLcLParticle);
      instance.SetDelete(&delete_physcLcLParticle);
      instance.SetDeleteArray(&deleteArray_physcLcLParticle);
      instance.SetDestructor(&destruct_physcLcLParticle);
      return &instance;
   }
   TGenericClassInfo *GenerateInitInstance(const ::phys::Particle*)
   {
      return GenerateInitInstanceLocal((::phys::Particle*)0);
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal((const ::phys::Particle*)0x0); R__UseDummy(_R__UNIQUE_DICT_(Init));
} // end of namespace ROOT

namespace ROOT {
   static void *new_physcLcLJet(void *p = 0);
   static void *newArray_physcLcLJet(Long_t size, void *p);
   static void delete_physcLcLJet(void *p);
   static void deleteArray_physcLcLJet(void *p);
   static void destruct_physcLcLJet(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const ::phys::Jet*)
   {
      ::phys::Jet *ptr = 0;
      static ::TVirtualIsAProxy* isa_proxy = new ::TInstrumentedIsAProxy< ::phys::Jet >(0);
      static ::ROOT::TGenericClassInfo 
         instance("phys::Jet", ::phys::Jet::Class_Version(), "Jet.h", 18,
                  typeid(::phys::Jet), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &::phys::Jet::Dictionary, isa_proxy, 4,
                  sizeof(::phys::Jet) );
      instance.SetNew(&new_physcLcLJet);
      instance.SetNewArray(&newArray_physcLcLJet);
      instance.SetDelete(&delete_physcLcLJet);
      instance.SetDeleteArray(&deleteArray_physcLcLJet);
      instance.SetDestructor(&destruct_physcLcLJet);
      return &instance;
   }
   TGenericClassInfo *GenerateInitInstance(const ::phys::Jet*)
   {
      return GenerateInitInstanceLocal((::phys::Jet*)0);
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal((const ::phys::Jet*)0x0); R__UseDummy(_R__UNIQUE_DICT_(Init));
} // end of namespace ROOT

namespace ROOT {
   static void *new_physcLcLLepton(void *p = 0);
   static void *newArray_physcLcLLepton(Long_t size, void *p);
   static void delete_physcLcLLepton(void *p);
   static void deleteArray_physcLcLLepton(void *p);
   static void destruct_physcLcLLepton(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const ::phys::Lepton*)
   {
      ::phys::Lepton *ptr = 0;
      static ::TVirtualIsAProxy* isa_proxy = new ::TInstrumentedIsAProxy< ::phys::Lepton >(0);
      static ::ROOT::TGenericClassInfo 
         instance("phys::Lepton", ::phys::Lepton::Class_Version(), "Lepton.h", 19,
                  typeid(::phys::Lepton), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &::phys::Lepton::Dictionary, isa_proxy, 4,
                  sizeof(::phys::Lepton) );
      instance.SetNew(&new_physcLcLLepton);
      instance.SetNewArray(&newArray_physcLcLLepton);
      instance.SetDelete(&delete_physcLcLLepton);
      instance.SetDeleteArray(&deleteArray_physcLcLLepton);
      instance.SetDestructor(&destruct_physcLcLLepton);
      return &instance;
   }
   TGenericClassInfo *GenerateInitInstance(const ::phys::Lepton*)
   {
      return GenerateInitInstanceLocal((::phys::Lepton*)0);
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal((const ::phys::Lepton*)0x0); R__UseDummy(_R__UNIQUE_DICT_(Init));
} // end of namespace ROOT

namespace ROOT {
   static void *new_physcLcLElectron(void *p = 0);
   static void *newArray_physcLcLElectron(Long_t size, void *p);
   static void delete_physcLcLElectron(void *p);
   static void deleteArray_physcLcLElectron(void *p);
   static void destruct_physcLcLElectron(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const ::phys::Electron*)
   {
      ::phys::Electron *ptr = 0;
      static ::TVirtualIsAProxy* isa_proxy = new ::TInstrumentedIsAProxy< ::phys::Electron >(0);
      static ::ROOT::TGenericClassInfo 
         instance("phys::Electron", ::phys::Electron::Class_Version(), "Electron.h", 16,
                  typeid(::phys::Electron), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &::phys::Electron::Dictionary, isa_proxy, 4,
                  sizeof(::phys::Electron) );
      instance.SetNew(&new_physcLcLElectron);
      instance.SetNewArray(&newArray_physcLcLElectron);
      instance.SetDelete(&delete_physcLcLElectron);
      instance.SetDeleteArray(&deleteArray_physcLcLElectron);
      instance.SetDestructor(&destruct_physcLcLElectron);
      return &instance;
   }
   TGenericClassInfo *GenerateInitInstance(const ::phys::Electron*)
   {
      return GenerateInitInstanceLocal((::phys::Electron*)0);
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal((const ::phys::Electron*)0x0); R__UseDummy(_R__UNIQUE_DICT_(Init));
} // end of namespace ROOT

namespace ROOT {
   static TClass *physcLcLBosonlEphyscLcLParticlegR_Dictionary();
   static void physcLcLBosonlEphyscLcLParticlegR_TClassManip(TClass*);
   static void *new_physcLcLBosonlEphyscLcLParticlegR(void *p = 0);
   static void *newArray_physcLcLBosonlEphyscLcLParticlegR(Long_t size, void *p);
   static void delete_physcLcLBosonlEphyscLcLParticlegR(void *p);
   static void deleteArray_physcLcLBosonlEphyscLcLParticlegR(void *p);
   static void destruct_physcLcLBosonlEphyscLcLParticlegR(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const ::phys::Boson<phys::Particle>*)
   {
      ::phys::Boson<phys::Particle> *ptr = 0;
      static ::TVirtualIsAProxy* isa_proxy = new ::TInstrumentedIsAProxy< ::phys::Boson<phys::Particle> >(0);
      static ::ROOT::TGenericClassInfo 
         instance("phys::Boson<phys::Particle>", ::phys::Boson<phys::Particle>::Class_Version(), "Boson.h", 18,
                  typeid(::phys::Boson<phys::Particle>), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &physcLcLBosonlEphyscLcLParticlegR_Dictionary, isa_proxy, 4,
                  sizeof(::phys::Boson<phys::Particle>) );
      instance.SetNew(&new_physcLcLBosonlEphyscLcLParticlegR);
      instance.SetNewArray(&newArray_physcLcLBosonlEphyscLcLParticlegR);
      instance.SetDelete(&delete_physcLcLBosonlEphyscLcLParticlegR);
      instance.SetDeleteArray(&deleteArray_physcLcLBosonlEphyscLcLParticlegR);
      instance.SetDestructor(&destruct_physcLcLBosonlEphyscLcLParticlegR);
      return &instance;
   }
   TGenericClassInfo *GenerateInitInstance(const ::phys::Boson<phys::Particle>*)
   {
      return GenerateInitInstanceLocal((::phys::Boson<phys::Particle>*)0);
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal((const ::phys::Boson<phys::Particle>*)0x0); R__UseDummy(_R__UNIQUE_DICT_(Init));

   // Dictionary for non-ClassDef classes
   static TClass *physcLcLBosonlEphyscLcLParticlegR_Dictionary() {
      TClass* theClass =::ROOT::GenerateInitInstanceLocal((const ::phys::Boson<phys::Particle>*)0x0)->GetClass();
      physcLcLBosonlEphyscLcLParticlegR_TClassManip(theClass);
   return theClass;
   }

   static void physcLcLBosonlEphyscLcLParticlegR_TClassManip(TClass* ){
   }

} // end of namespace ROOT

namespace ROOT {
   static TClass *physcLcLBosonlEphyscLcLLeptongR_Dictionary();
   static void physcLcLBosonlEphyscLcLLeptongR_TClassManip(TClass*);
   static void *new_physcLcLBosonlEphyscLcLLeptongR(void *p = 0);
   static void *newArray_physcLcLBosonlEphyscLcLLeptongR(Long_t size, void *p);
   static void delete_physcLcLBosonlEphyscLcLLeptongR(void *p);
   static void deleteArray_physcLcLBosonlEphyscLcLLeptongR(void *p);
   static void destruct_physcLcLBosonlEphyscLcLLeptongR(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const ::phys::Boson<phys::Lepton>*)
   {
      ::phys::Boson<phys::Lepton> *ptr = 0;
      static ::TVirtualIsAProxy* isa_proxy = new ::TInstrumentedIsAProxy< ::phys::Boson<phys::Lepton> >(0);
      static ::ROOT::TGenericClassInfo 
         instance("phys::Boson<phys::Lepton>", ::phys::Boson<phys::Lepton>::Class_Version(), "Boson.h", 18,
                  typeid(::phys::Boson<phys::Lepton>), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &physcLcLBosonlEphyscLcLLeptongR_Dictionary, isa_proxy, 4,
                  sizeof(::phys::Boson<phys::Lepton>) );
      instance.SetNew(&new_physcLcLBosonlEphyscLcLLeptongR);
      instance.SetNewArray(&newArray_physcLcLBosonlEphyscLcLLeptongR);
      instance.SetDelete(&delete_physcLcLBosonlEphyscLcLLeptongR);
      instance.SetDeleteArray(&deleteArray_physcLcLBosonlEphyscLcLLeptongR);
      instance.SetDestructor(&destruct_physcLcLBosonlEphyscLcLLeptongR);
      return &instance;
   }
   TGenericClassInfo *GenerateInitInstance(const ::phys::Boson<phys::Lepton>*)
   {
      return GenerateInitInstanceLocal((::phys::Boson<phys::Lepton>*)0);
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal((const ::phys::Boson<phys::Lepton>*)0x0); R__UseDummy(_R__UNIQUE_DICT_(Init));

   // Dictionary for non-ClassDef classes
   static TClass *physcLcLBosonlEphyscLcLLeptongR_Dictionary() {
      TClass* theClass =::ROOT::GenerateInitInstanceLocal((const ::phys::Boson<phys::Lepton>*)0x0)->GetClass();
      physcLcLBosonlEphyscLcLLeptongR_TClassManip(theClass);
   return theClass;
   }

   static void physcLcLBosonlEphyscLcLLeptongR_TClassManip(TClass* ){
   }

} // end of namespace ROOT

namespace ROOT {
   static TClass *physcLcLBosonlEphyscLcLJetgR_Dictionary();
   static void physcLcLBosonlEphyscLcLJetgR_TClassManip(TClass*);
   static void *new_physcLcLBosonlEphyscLcLJetgR(void *p = 0);
   static void *newArray_physcLcLBosonlEphyscLcLJetgR(Long_t size, void *p);
   static void delete_physcLcLBosonlEphyscLcLJetgR(void *p);
   static void deleteArray_physcLcLBosonlEphyscLcLJetgR(void *p);
   static void destruct_physcLcLBosonlEphyscLcLJetgR(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const ::phys::Boson<phys::Jet>*)
   {
      ::phys::Boson<phys::Jet> *ptr = 0;
      static ::TVirtualIsAProxy* isa_proxy = new ::TInstrumentedIsAProxy< ::phys::Boson<phys::Jet> >(0);
      static ::ROOT::TGenericClassInfo 
         instance("phys::Boson<phys::Jet>", ::phys::Boson<phys::Jet>::Class_Version(), "Boson.h", 18,
                  typeid(::phys::Boson<phys::Jet>), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &physcLcLBosonlEphyscLcLJetgR_Dictionary, isa_proxy, 4,
                  sizeof(::phys::Boson<phys::Jet>) );
      instance.SetNew(&new_physcLcLBosonlEphyscLcLJetgR);
      instance.SetNewArray(&newArray_physcLcLBosonlEphyscLcLJetgR);
      instance.SetDelete(&delete_physcLcLBosonlEphyscLcLJetgR);
      instance.SetDeleteArray(&deleteArray_physcLcLBosonlEphyscLcLJetgR);
      instance.SetDestructor(&destruct_physcLcLBosonlEphyscLcLJetgR);
      return &instance;
   }
   TGenericClassInfo *GenerateInitInstance(const ::phys::Boson<phys::Jet>*)
   {
      return GenerateInitInstanceLocal((::phys::Boson<phys::Jet>*)0);
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal((const ::phys::Boson<phys::Jet>*)0x0); R__UseDummy(_R__UNIQUE_DICT_(Init));

   // Dictionary for non-ClassDef classes
   static TClass *physcLcLBosonlEphyscLcLJetgR_Dictionary() {
      TClass* theClass =::ROOT::GenerateInitInstanceLocal((const ::phys::Boson<phys::Jet>*)0x0)->GetClass();
      physcLcLBosonlEphyscLcLJetgR_TClassManip(theClass);
   return theClass;
   }

   static void physcLcLBosonlEphyscLcLJetgR_TClassManip(TClass* ){
   }

} // end of namespace ROOT

namespace ROOT {
   static TClass *physcLcLBosonlEphyscLcLElectrongR_Dictionary();
   static void physcLcLBosonlEphyscLcLElectrongR_TClassManip(TClass*);
   static void *new_physcLcLBosonlEphyscLcLElectrongR(void *p = 0);
   static void *newArray_physcLcLBosonlEphyscLcLElectrongR(Long_t size, void *p);
   static void delete_physcLcLBosonlEphyscLcLElectrongR(void *p);
   static void deleteArray_physcLcLBosonlEphyscLcLElectrongR(void *p);
   static void destruct_physcLcLBosonlEphyscLcLElectrongR(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const ::phys::Boson<phys::Electron>*)
   {
      ::phys::Boson<phys::Electron> *ptr = 0;
      static ::TVirtualIsAProxy* isa_proxy = new ::TInstrumentedIsAProxy< ::phys::Boson<phys::Electron> >(0);
      static ::ROOT::TGenericClassInfo 
         instance("phys::Boson<phys::Electron>", ::phys::Boson<phys::Electron>::Class_Version(), "Boson.h", 18,
                  typeid(::phys::Boson<phys::Electron>), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &physcLcLBosonlEphyscLcLElectrongR_Dictionary, isa_proxy, 4,
                  sizeof(::phys::Boson<phys::Electron>) );
      instance.SetNew(&new_physcLcLBosonlEphyscLcLElectrongR);
      instance.SetNewArray(&newArray_physcLcLBosonlEphyscLcLElectrongR);
      instance.SetDelete(&delete_physcLcLBosonlEphyscLcLElectrongR);
      instance.SetDeleteArray(&deleteArray_physcLcLBosonlEphyscLcLElectrongR);
      instance.SetDestructor(&destruct_physcLcLBosonlEphyscLcLElectrongR);
      return &instance;
   }
   TGenericClassInfo *GenerateInitInstance(const ::phys::Boson<phys::Electron>*)
   {
      return GenerateInitInstanceLocal((::phys::Boson<phys::Electron>*)0);
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal((const ::phys::Boson<phys::Electron>*)0x0); R__UseDummy(_R__UNIQUE_DICT_(Init));

   // Dictionary for non-ClassDef classes
   static TClass *physcLcLBosonlEphyscLcLElectrongR_Dictionary() {
      TClass* theClass =::ROOT::GenerateInitInstanceLocal((const ::phys::Boson<phys::Electron>*)0x0)->GetClass();
      physcLcLBosonlEphyscLcLElectrongR_TClassManip(theClass);
   return theClass;
   }

   static void physcLcLBosonlEphyscLcLElectrongR_TClassManip(TClass* ){
   }

} // end of namespace ROOT

namespace ROOT {
   static TClass *physcLcLDiBosonlEphyscLcLLeptoncOphyscLcLLeptongR_Dictionary();
   static void physcLcLDiBosonlEphyscLcLLeptoncOphyscLcLLeptongR_TClassManip(TClass*);
   static void *new_physcLcLDiBosonlEphyscLcLLeptoncOphyscLcLLeptongR(void *p = 0);
   static void *newArray_physcLcLDiBosonlEphyscLcLLeptoncOphyscLcLLeptongR(Long_t size, void *p);
   static void delete_physcLcLDiBosonlEphyscLcLLeptoncOphyscLcLLeptongR(void *p);
   static void deleteArray_physcLcLDiBosonlEphyscLcLLeptoncOphyscLcLLeptongR(void *p);
   static void destruct_physcLcLDiBosonlEphyscLcLLeptoncOphyscLcLLeptongR(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const ::phys::DiBoson<phys::Lepton,phys::Lepton>*)
   {
      ::phys::DiBoson<phys::Lepton,phys::Lepton> *ptr = 0;
      static ::TVirtualIsAProxy* isa_proxy = new ::TInstrumentedIsAProxy< ::phys::DiBoson<phys::Lepton,phys::Lepton> >(0);
      static ::ROOT::TGenericClassInfo 
         instance("phys::DiBoson<phys::Lepton,phys::Lepton>", ::phys::DiBoson<phys::Lepton,phys::Lepton>::Class_Version(), "DiBoson.h", 19,
                  typeid(::phys::DiBoson<phys::Lepton,phys::Lepton>), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &physcLcLDiBosonlEphyscLcLLeptoncOphyscLcLLeptongR_Dictionary, isa_proxy, 4,
                  sizeof(::phys::DiBoson<phys::Lepton,phys::Lepton>) );
      instance.SetNew(&new_physcLcLDiBosonlEphyscLcLLeptoncOphyscLcLLeptongR);
      instance.SetNewArray(&newArray_physcLcLDiBosonlEphyscLcLLeptoncOphyscLcLLeptongR);
      instance.SetDelete(&delete_physcLcLDiBosonlEphyscLcLLeptoncOphyscLcLLeptongR);
      instance.SetDeleteArray(&deleteArray_physcLcLDiBosonlEphyscLcLLeptoncOphyscLcLLeptongR);
      instance.SetDestructor(&destruct_physcLcLDiBosonlEphyscLcLLeptoncOphyscLcLLeptongR);
      return &instance;
   }
   TGenericClassInfo *GenerateInitInstance(const ::phys::DiBoson<phys::Lepton,phys::Lepton>*)
   {
      return GenerateInitInstanceLocal((::phys::DiBoson<phys::Lepton,phys::Lepton>*)0);
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal((const ::phys::DiBoson<phys::Lepton,phys::Lepton>*)0x0); R__UseDummy(_R__UNIQUE_DICT_(Init));

   // Dictionary for non-ClassDef classes
   static TClass *physcLcLDiBosonlEphyscLcLLeptoncOphyscLcLLeptongR_Dictionary() {
      TClass* theClass =::ROOT::GenerateInitInstanceLocal((const ::phys::DiBoson<phys::Lepton,phys::Lepton>*)0x0)->GetClass();
      physcLcLDiBosonlEphyscLcLLeptoncOphyscLcLLeptongR_TClassManip(theClass);
   return theClass;
   }

   static void physcLcLDiBosonlEphyscLcLLeptoncOphyscLcLLeptongR_TClassManip(TClass* ){
   }

} // end of namespace ROOT

namespace ROOT {
   static TClass *physcLcLDiBosonlEphyscLcLElectroncOphyscLcLElectrongR_Dictionary();
   static void physcLcLDiBosonlEphyscLcLElectroncOphyscLcLElectrongR_TClassManip(TClass*);
   static void *new_physcLcLDiBosonlEphyscLcLElectroncOphyscLcLElectrongR(void *p = 0);
   static void *newArray_physcLcLDiBosonlEphyscLcLElectroncOphyscLcLElectrongR(Long_t size, void *p);
   static void delete_physcLcLDiBosonlEphyscLcLElectroncOphyscLcLElectrongR(void *p);
   static void deleteArray_physcLcLDiBosonlEphyscLcLElectroncOphyscLcLElectrongR(void *p);
   static void destruct_physcLcLDiBosonlEphyscLcLElectroncOphyscLcLElectrongR(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const ::phys::DiBoson<phys::Electron,phys::Electron>*)
   {
      ::phys::DiBoson<phys::Electron,phys::Electron> *ptr = 0;
      static ::TVirtualIsAProxy* isa_proxy = new ::TInstrumentedIsAProxy< ::phys::DiBoson<phys::Electron,phys::Electron> >(0);
      static ::ROOT::TGenericClassInfo 
         instance("phys::DiBoson<phys::Electron,phys::Electron>", ::phys::DiBoson<phys::Electron,phys::Electron>::Class_Version(), "DiBoson.h", 19,
                  typeid(::phys::DiBoson<phys::Electron,phys::Electron>), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &physcLcLDiBosonlEphyscLcLElectroncOphyscLcLElectrongR_Dictionary, isa_proxy, 4,
                  sizeof(::phys::DiBoson<phys::Electron,phys::Electron>) );
      instance.SetNew(&new_physcLcLDiBosonlEphyscLcLElectroncOphyscLcLElectrongR);
      instance.SetNewArray(&newArray_physcLcLDiBosonlEphyscLcLElectroncOphyscLcLElectrongR);
      instance.SetDelete(&delete_physcLcLDiBosonlEphyscLcLElectroncOphyscLcLElectrongR);
      instance.SetDeleteArray(&deleteArray_physcLcLDiBosonlEphyscLcLElectroncOphyscLcLElectrongR);
      instance.SetDestructor(&destruct_physcLcLDiBosonlEphyscLcLElectroncOphyscLcLElectrongR);
      return &instance;
   }
   TGenericClassInfo *GenerateInitInstance(const ::phys::DiBoson<phys::Electron,phys::Electron>*)
   {
      return GenerateInitInstanceLocal((::phys::DiBoson<phys::Electron,phys::Electron>*)0);
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal((const ::phys::DiBoson<phys::Electron,phys::Electron>*)0x0); R__UseDummy(_R__UNIQUE_DICT_(Init));

   // Dictionary for non-ClassDef classes
   static TClass *physcLcLDiBosonlEphyscLcLElectroncOphyscLcLElectrongR_Dictionary() {
      TClass* theClass =::ROOT::GenerateInitInstanceLocal((const ::phys::DiBoson<phys::Electron,phys::Electron>*)0x0)->GetClass();
      physcLcLDiBosonlEphyscLcLElectroncOphyscLcLElectrongR_TClassManip(theClass);
   return theClass;
   }

   static void physcLcLDiBosonlEphyscLcLElectroncOphyscLcLElectrongR_TClassManip(TClass* ){
   }

} // end of namespace ROOT

namespace ROOT {
   static TClass *physcLcLDiBosonlEphyscLcLElectroncOphyscLcLLeptongR_Dictionary();
   static void physcLcLDiBosonlEphyscLcLElectroncOphyscLcLLeptongR_TClassManip(TClass*);
   static void *new_physcLcLDiBosonlEphyscLcLElectroncOphyscLcLLeptongR(void *p = 0);
   static void *newArray_physcLcLDiBosonlEphyscLcLElectroncOphyscLcLLeptongR(Long_t size, void *p);
   static void delete_physcLcLDiBosonlEphyscLcLElectroncOphyscLcLLeptongR(void *p);
   static void deleteArray_physcLcLDiBosonlEphyscLcLElectroncOphyscLcLLeptongR(void *p);
   static void destruct_physcLcLDiBosonlEphyscLcLElectroncOphyscLcLLeptongR(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const ::phys::DiBoson<phys::Electron,phys::Lepton>*)
   {
      ::phys::DiBoson<phys::Electron,phys::Lepton> *ptr = 0;
      static ::TVirtualIsAProxy* isa_proxy = new ::TInstrumentedIsAProxy< ::phys::DiBoson<phys::Electron,phys::Lepton> >(0);
      static ::ROOT::TGenericClassInfo 
         instance("phys::DiBoson<phys::Electron,phys::Lepton>", ::phys::DiBoson<phys::Electron,phys::Lepton>::Class_Version(), "DiBoson.h", 19,
                  typeid(::phys::DiBoson<phys::Electron,phys::Lepton>), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &physcLcLDiBosonlEphyscLcLElectroncOphyscLcLLeptongR_Dictionary, isa_proxy, 4,
                  sizeof(::phys::DiBoson<phys::Electron,phys::Lepton>) );
      instance.SetNew(&new_physcLcLDiBosonlEphyscLcLElectroncOphyscLcLLeptongR);
      instance.SetNewArray(&newArray_physcLcLDiBosonlEphyscLcLElectroncOphyscLcLLeptongR);
      instance.SetDelete(&delete_physcLcLDiBosonlEphyscLcLElectroncOphyscLcLLeptongR);
      instance.SetDeleteArray(&deleteArray_physcLcLDiBosonlEphyscLcLElectroncOphyscLcLLeptongR);
      instance.SetDestructor(&destruct_physcLcLDiBosonlEphyscLcLElectroncOphyscLcLLeptongR);
      return &instance;
   }
   TGenericClassInfo *GenerateInitInstance(const ::phys::DiBoson<phys::Electron,phys::Lepton>*)
   {
      return GenerateInitInstanceLocal((::phys::DiBoson<phys::Electron,phys::Lepton>*)0);
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal((const ::phys::DiBoson<phys::Electron,phys::Lepton>*)0x0); R__UseDummy(_R__UNIQUE_DICT_(Init));

   // Dictionary for non-ClassDef classes
   static TClass *physcLcLDiBosonlEphyscLcLElectroncOphyscLcLLeptongR_Dictionary() {
      TClass* theClass =::ROOT::GenerateInitInstanceLocal((const ::phys::DiBoson<phys::Electron,phys::Lepton>*)0x0)->GetClass();
      physcLcLDiBosonlEphyscLcLElectroncOphyscLcLLeptongR_TClassManip(theClass);
   return theClass;
   }

   static void physcLcLDiBosonlEphyscLcLElectroncOphyscLcLLeptongR_TClassManip(TClass* ){
   }

} // end of namespace ROOT

namespace ROOT {
   static TClass *physcLcLDiBosonlEphyscLcLLeptoncOphyscLcLElectrongR_Dictionary();
   static void physcLcLDiBosonlEphyscLcLLeptoncOphyscLcLElectrongR_TClassManip(TClass*);
   static void *new_physcLcLDiBosonlEphyscLcLLeptoncOphyscLcLElectrongR(void *p = 0);
   static void *newArray_physcLcLDiBosonlEphyscLcLLeptoncOphyscLcLElectrongR(Long_t size, void *p);
   static void delete_physcLcLDiBosonlEphyscLcLLeptoncOphyscLcLElectrongR(void *p);
   static void deleteArray_physcLcLDiBosonlEphyscLcLLeptoncOphyscLcLElectrongR(void *p);
   static void destruct_physcLcLDiBosonlEphyscLcLLeptoncOphyscLcLElectrongR(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const ::phys::DiBoson<phys::Lepton,phys::Electron>*)
   {
      ::phys::DiBoson<phys::Lepton,phys::Electron> *ptr = 0;
      static ::TVirtualIsAProxy* isa_proxy = new ::TInstrumentedIsAProxy< ::phys::DiBoson<phys::Lepton,phys::Electron> >(0);
      static ::ROOT::TGenericClassInfo 
         instance("phys::DiBoson<phys::Lepton,phys::Electron>", ::phys::DiBoson<phys::Lepton,phys::Electron>::Class_Version(), "DiBoson.h", 19,
                  typeid(::phys::DiBoson<phys::Lepton,phys::Electron>), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &physcLcLDiBosonlEphyscLcLLeptoncOphyscLcLElectrongR_Dictionary, isa_proxy, 4,
                  sizeof(::phys::DiBoson<phys::Lepton,phys::Electron>) );
      instance.SetNew(&new_physcLcLDiBosonlEphyscLcLLeptoncOphyscLcLElectrongR);
      instance.SetNewArray(&newArray_physcLcLDiBosonlEphyscLcLLeptoncOphyscLcLElectrongR);
      instance.SetDelete(&delete_physcLcLDiBosonlEphyscLcLLeptoncOphyscLcLElectrongR);
      instance.SetDeleteArray(&deleteArray_physcLcLDiBosonlEphyscLcLLeptoncOphyscLcLElectrongR);
      instance.SetDestructor(&destruct_physcLcLDiBosonlEphyscLcLLeptoncOphyscLcLElectrongR);
      return &instance;
   }
   TGenericClassInfo *GenerateInitInstance(const ::phys::DiBoson<phys::Lepton,phys::Electron>*)
   {
      return GenerateInitInstanceLocal((::phys::DiBoson<phys::Lepton,phys::Electron>*)0);
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal((const ::phys::DiBoson<phys::Lepton,phys::Electron>*)0x0); R__UseDummy(_R__UNIQUE_DICT_(Init));

   // Dictionary for non-ClassDef classes
   static TClass *physcLcLDiBosonlEphyscLcLLeptoncOphyscLcLElectrongR_Dictionary() {
      TClass* theClass =::ROOT::GenerateInitInstanceLocal((const ::phys::DiBoson<phys::Lepton,phys::Electron>*)0x0)->GetClass();
      physcLcLDiBosonlEphyscLcLLeptoncOphyscLcLElectrongR_TClassManip(theClass);
   return theClass;
   }

   static void physcLcLDiBosonlEphyscLcLLeptoncOphyscLcLElectrongR_TClassManip(TClass* ){
   }

} // end of namespace ROOT

namespace ROOT {
   static TClass *physcLcLDiBosonlEphyscLcLParticlecOphyscLcLParticlegR_Dictionary();
   static void physcLcLDiBosonlEphyscLcLParticlecOphyscLcLParticlegR_TClassManip(TClass*);
   static void *new_physcLcLDiBosonlEphyscLcLParticlecOphyscLcLParticlegR(void *p = 0);
   static void *newArray_physcLcLDiBosonlEphyscLcLParticlecOphyscLcLParticlegR(Long_t size, void *p);
   static void delete_physcLcLDiBosonlEphyscLcLParticlecOphyscLcLParticlegR(void *p);
   static void deleteArray_physcLcLDiBosonlEphyscLcLParticlecOphyscLcLParticlegR(void *p);
   static void destruct_physcLcLDiBosonlEphyscLcLParticlecOphyscLcLParticlegR(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const ::phys::DiBoson<phys::Particle,phys::Particle>*)
   {
      ::phys::DiBoson<phys::Particle,phys::Particle> *ptr = 0;
      static ::TVirtualIsAProxy* isa_proxy = new ::TInstrumentedIsAProxy< ::phys::DiBoson<phys::Particle,phys::Particle> >(0);
      static ::ROOT::TGenericClassInfo 
         instance("phys::DiBoson<phys::Particle,phys::Particle>", ::phys::DiBoson<phys::Particle,phys::Particle>::Class_Version(), "DiBoson.h", 19,
                  typeid(::phys::DiBoson<phys::Particle,phys::Particle>), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &physcLcLDiBosonlEphyscLcLParticlecOphyscLcLParticlegR_Dictionary, isa_proxy, 4,
                  sizeof(::phys::DiBoson<phys::Particle,phys::Particle>) );
      instance.SetNew(&new_physcLcLDiBosonlEphyscLcLParticlecOphyscLcLParticlegR);
      instance.SetNewArray(&newArray_physcLcLDiBosonlEphyscLcLParticlecOphyscLcLParticlegR);
      instance.SetDelete(&delete_physcLcLDiBosonlEphyscLcLParticlecOphyscLcLParticlegR);
      instance.SetDeleteArray(&deleteArray_physcLcLDiBosonlEphyscLcLParticlecOphyscLcLParticlegR);
      instance.SetDestructor(&destruct_physcLcLDiBosonlEphyscLcLParticlecOphyscLcLParticlegR);
      return &instance;
   }
   TGenericClassInfo *GenerateInitInstance(const ::phys::DiBoson<phys::Particle,phys::Particle>*)
   {
      return GenerateInitInstanceLocal((::phys::DiBoson<phys::Particle,phys::Particle>*)0);
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal((const ::phys::DiBoson<phys::Particle,phys::Particle>*)0x0); R__UseDummy(_R__UNIQUE_DICT_(Init));

   // Dictionary for non-ClassDef classes
   static TClass *physcLcLDiBosonlEphyscLcLParticlecOphyscLcLParticlegR_Dictionary() {
      TClass* theClass =::ROOT::GenerateInitInstanceLocal((const ::phys::DiBoson<phys::Particle,phys::Particle>*)0x0)->GetClass();
      physcLcLDiBosonlEphyscLcLParticlecOphyscLcLParticlegR_TClassManip(theClass);
   return theClass;
   }

   static void physcLcLDiBosonlEphyscLcLParticlecOphyscLcLParticlegR_TClassManip(TClass* ){
   }

} // end of namespace ROOT

namespace ROOT {
   static void *new_physcLcLGenEventWeights(void *p = 0);
   static void *newArray_physcLcLGenEventWeights(Long_t size, void *p);
   static void delete_physcLcLGenEventWeights(void *p);
   static void deleteArray_physcLcLGenEventWeights(void *p);
   static void destruct_physcLcLGenEventWeights(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const ::phys::GenEventWeights*)
   {
      ::phys::GenEventWeights *ptr = 0;
      static ::TVirtualIsAProxy* isa_proxy = new ::TInstrumentedIsAProxy< ::phys::GenEventWeights >(0);
      static ::ROOT::TGenericClassInfo 
         instance("phys::GenEventWeights", ::phys::GenEventWeights::Class_Version(), "GenEventWeights.h", 22,
                  typeid(::phys::GenEventWeights), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &::phys::GenEventWeights::Dictionary, isa_proxy, 4,
                  sizeof(::phys::GenEventWeights) );
      instance.SetNew(&new_physcLcLGenEventWeights);
      instance.SetNewArray(&newArray_physcLcLGenEventWeights);
      instance.SetDelete(&delete_physcLcLGenEventWeights);
      instance.SetDeleteArray(&deleteArray_physcLcLGenEventWeights);
      instance.SetDestructor(&destruct_physcLcLGenEventWeights);
      return &instance;
   }
   TGenericClassInfo *GenerateInitInstance(const ::phys::GenEventWeights*)
   {
      return GenerateInitInstanceLocal((::phys::GenEventWeights*)0);
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal((const ::phys::GenEventWeights*)0x0); R__UseDummy(_R__UNIQUE_DICT_(Init));
} // end of namespace ROOT

namespace ROOT {
   static void *new_physcLcLMELA(void *p = 0);
   static void *newArray_physcLcLMELA(Long_t size, void *p);
   static void delete_physcLcLMELA(void *p);
   static void deleteArray_physcLcLMELA(void *p);
   static void destruct_physcLcLMELA(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const ::phys::MELA*)
   {
      ::phys::MELA *ptr = 0;
      static ::TVirtualIsAProxy* isa_proxy = new ::TInstrumentedIsAProxy< ::phys::MELA >(0);
      static ::ROOT::TGenericClassInfo 
         instance("phys::MELA", ::phys::MELA::Class_Version(), "MELA.h", 22,
                  typeid(::phys::MELA), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &::phys::MELA::Dictionary, isa_proxy, 4,
                  sizeof(::phys::MELA) );
      instance.SetNew(&new_physcLcLMELA);
      instance.SetNewArray(&newArray_physcLcLMELA);
      instance.SetDelete(&delete_physcLcLMELA);
      instance.SetDeleteArray(&deleteArray_physcLcLMELA);
      instance.SetDestructor(&destruct_physcLcLMELA);
      return &instance;
   }
   TGenericClassInfo *GenerateInitInstance(const ::phys::MELA*)
   {
      return GenerateInitInstanceLocal((::phys::MELA*)0);
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal((const ::phys::MELA*)0x0); R__UseDummy(_R__UNIQUE_DICT_(Init));
} // end of namespace ROOT

namespace phys {
//______________________________________________________________________________
atomic_TClass_ptr Particle::fgIsA(0);  // static to hold class pointer

//______________________________________________________________________________
const char *Particle::Class_Name()
{
   return "phys::Particle";
}

//______________________________________________________________________________
const char *Particle::ImplFileName()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::phys::Particle*)0x0)->GetImplFileName();
}

//______________________________________________________________________________
int Particle::ImplFileLine()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::phys::Particle*)0x0)->GetImplFileLine();
}

//______________________________________________________________________________
TClass *Particle::Dictionary()
{
   fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::phys::Particle*)0x0)->GetClass();
   return fgIsA;
}

//______________________________________________________________________________
TClass *Particle::Class()
{
   if (!fgIsA.load()) { R__LOCKGUARD(gInterpreterMutex); fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::phys::Particle*)0x0)->GetClass(); }
   return fgIsA;
}

} // namespace phys
namespace phys {
//______________________________________________________________________________
atomic_TClass_ptr Jet::fgIsA(0);  // static to hold class pointer

//______________________________________________________________________________
const char *Jet::Class_Name()
{
   return "phys::Jet";
}

//______________________________________________________________________________
const char *Jet::ImplFileName()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::phys::Jet*)0x0)->GetImplFileName();
}

//______________________________________________________________________________
int Jet::ImplFileLine()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::phys::Jet*)0x0)->GetImplFileLine();
}

//______________________________________________________________________________
TClass *Jet::Dictionary()
{
   fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::phys::Jet*)0x0)->GetClass();
   return fgIsA;
}

//______________________________________________________________________________
TClass *Jet::Class()
{
   if (!fgIsA.load()) { R__LOCKGUARD(gInterpreterMutex); fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::phys::Jet*)0x0)->GetClass(); }
   return fgIsA;
}

} // namespace phys
namespace phys {
//______________________________________________________________________________
atomic_TClass_ptr Lepton::fgIsA(0);  // static to hold class pointer

//______________________________________________________________________________
const char *Lepton::Class_Name()
{
   return "phys::Lepton";
}

//______________________________________________________________________________
const char *Lepton::ImplFileName()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::phys::Lepton*)0x0)->GetImplFileName();
}

//______________________________________________________________________________
int Lepton::ImplFileLine()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::phys::Lepton*)0x0)->GetImplFileLine();
}

//______________________________________________________________________________
TClass *Lepton::Dictionary()
{
   fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::phys::Lepton*)0x0)->GetClass();
   return fgIsA;
}

//______________________________________________________________________________
TClass *Lepton::Class()
{
   if (!fgIsA.load()) { R__LOCKGUARD(gInterpreterMutex); fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::phys::Lepton*)0x0)->GetClass(); }
   return fgIsA;
}

} // namespace phys
namespace phys {
//______________________________________________________________________________
atomic_TClass_ptr Electron::fgIsA(0);  // static to hold class pointer

//______________________________________________________________________________
const char *Electron::Class_Name()
{
   return "phys::Electron";
}

//______________________________________________________________________________
const char *Electron::ImplFileName()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::phys::Electron*)0x0)->GetImplFileName();
}

//______________________________________________________________________________
int Electron::ImplFileLine()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::phys::Electron*)0x0)->GetImplFileLine();
}

//______________________________________________________________________________
TClass *Electron::Dictionary()
{
   fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::phys::Electron*)0x0)->GetClass();
   return fgIsA;
}

//______________________________________________________________________________
TClass *Electron::Class()
{
   if (!fgIsA.load()) { R__LOCKGUARD(gInterpreterMutex); fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::phys::Electron*)0x0)->GetClass(); }
   return fgIsA;
}

} // namespace phys
namespace phys {
//______________________________________________________________________________
template <> atomic_TClass_ptr Boson<phys::Particle>::fgIsA(0);  // static to hold class pointer

//______________________________________________________________________________
template <> const char *Boson<phys::Particle>::Class_Name()
{
   return "phys::Boson<phys::Particle>";
}

//______________________________________________________________________________
template <> const char *Boson<phys::Particle>::ImplFileName()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::phys::Boson<phys::Particle>*)0x0)->GetImplFileName();
}

//______________________________________________________________________________
template <> int Boson<phys::Particle>::ImplFileLine()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::phys::Boson<phys::Particle>*)0x0)->GetImplFileLine();
}

//______________________________________________________________________________
template <> TClass *Boson<phys::Particle>::Dictionary()
{
   fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::phys::Boson<phys::Particle>*)0x0)->GetClass();
   return fgIsA;
}

//______________________________________________________________________________
template <> TClass *Boson<phys::Particle>::Class()
{
   if (!fgIsA.load()) { R__LOCKGUARD(gInterpreterMutex); fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::phys::Boson<phys::Particle>*)0x0)->GetClass(); }
   return fgIsA;
}

} // namespace phys
namespace phys {
//______________________________________________________________________________
template <> atomic_TClass_ptr Boson<phys::Lepton>::fgIsA(0);  // static to hold class pointer

//______________________________________________________________________________
template <> const char *Boson<phys::Lepton>::Class_Name()
{
   return "phys::Boson<phys::Lepton>";
}

//______________________________________________________________________________
template <> const char *Boson<phys::Lepton>::ImplFileName()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::phys::Boson<phys::Lepton>*)0x0)->GetImplFileName();
}

//______________________________________________________________________________
template <> int Boson<phys::Lepton>::ImplFileLine()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::phys::Boson<phys::Lepton>*)0x0)->GetImplFileLine();
}

//______________________________________________________________________________
template <> TClass *Boson<phys::Lepton>::Dictionary()
{
   fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::phys::Boson<phys::Lepton>*)0x0)->GetClass();
   return fgIsA;
}

//______________________________________________________________________________
template <> TClass *Boson<phys::Lepton>::Class()
{
   if (!fgIsA.load()) { R__LOCKGUARD(gInterpreterMutex); fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::phys::Boson<phys::Lepton>*)0x0)->GetClass(); }
   return fgIsA;
}

} // namespace phys
namespace phys {
//______________________________________________________________________________
template <> atomic_TClass_ptr Boson<phys::Jet>::fgIsA(0);  // static to hold class pointer

//______________________________________________________________________________
template <> const char *Boson<phys::Jet>::Class_Name()
{
   return "phys::Boson<phys::Jet>";
}

//______________________________________________________________________________
template <> const char *Boson<phys::Jet>::ImplFileName()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::phys::Boson<phys::Jet>*)0x0)->GetImplFileName();
}

//______________________________________________________________________________
template <> int Boson<phys::Jet>::ImplFileLine()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::phys::Boson<phys::Jet>*)0x0)->GetImplFileLine();
}

//______________________________________________________________________________
template <> TClass *Boson<phys::Jet>::Dictionary()
{
   fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::phys::Boson<phys::Jet>*)0x0)->GetClass();
   return fgIsA;
}

//______________________________________________________________________________
template <> TClass *Boson<phys::Jet>::Class()
{
   if (!fgIsA.load()) { R__LOCKGUARD(gInterpreterMutex); fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::phys::Boson<phys::Jet>*)0x0)->GetClass(); }
   return fgIsA;
}

} // namespace phys
namespace phys {
//______________________________________________________________________________
template <> atomic_TClass_ptr Boson<phys::Electron>::fgIsA(0);  // static to hold class pointer

//______________________________________________________________________________
template <> const char *Boson<phys::Electron>::Class_Name()
{
   return "phys::Boson<phys::Electron>";
}

//______________________________________________________________________________
template <> const char *Boson<phys::Electron>::ImplFileName()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::phys::Boson<phys::Electron>*)0x0)->GetImplFileName();
}

//______________________________________________________________________________
template <> int Boson<phys::Electron>::ImplFileLine()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::phys::Boson<phys::Electron>*)0x0)->GetImplFileLine();
}

//______________________________________________________________________________
template <> TClass *Boson<phys::Electron>::Dictionary()
{
   fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::phys::Boson<phys::Electron>*)0x0)->GetClass();
   return fgIsA;
}

//______________________________________________________________________________
template <> TClass *Boson<phys::Electron>::Class()
{
   if (!fgIsA.load()) { R__LOCKGUARD(gInterpreterMutex); fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::phys::Boson<phys::Electron>*)0x0)->GetClass(); }
   return fgIsA;
}

} // namespace phys
namespace phys {
//______________________________________________________________________________
template <> atomic_TClass_ptr DiBoson<phys::Lepton,phys::Lepton>::fgIsA(0);  // static to hold class pointer

//______________________________________________________________________________
template <> const char *DiBoson<phys::Lepton,phys::Lepton>::Class_Name()
{
   return "phys::DiBoson<phys::Lepton,phys::Lepton>";
}

//______________________________________________________________________________
template <> const char *DiBoson<phys::Lepton,phys::Lepton>::ImplFileName()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::phys::DiBoson<phys::Lepton,phys::Lepton>*)0x0)->GetImplFileName();
}

//______________________________________________________________________________
template <> int DiBoson<phys::Lepton,phys::Lepton>::ImplFileLine()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::phys::DiBoson<phys::Lepton,phys::Lepton>*)0x0)->GetImplFileLine();
}

//______________________________________________________________________________
template <> TClass *DiBoson<phys::Lepton,phys::Lepton>::Dictionary()
{
   fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::phys::DiBoson<phys::Lepton,phys::Lepton>*)0x0)->GetClass();
   return fgIsA;
}

//______________________________________________________________________________
template <> TClass *DiBoson<phys::Lepton,phys::Lepton>::Class()
{
   if (!fgIsA.load()) { R__LOCKGUARD(gInterpreterMutex); fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::phys::DiBoson<phys::Lepton,phys::Lepton>*)0x0)->GetClass(); }
   return fgIsA;
}

} // namespace phys
namespace phys {
//______________________________________________________________________________
template <> atomic_TClass_ptr DiBoson<phys::Electron,phys::Electron>::fgIsA(0);  // static to hold class pointer

//______________________________________________________________________________
template <> const char *DiBoson<phys::Electron,phys::Electron>::Class_Name()
{
   return "phys::DiBoson<phys::Electron,phys::Electron>";
}

//______________________________________________________________________________
template <> const char *DiBoson<phys::Electron,phys::Electron>::ImplFileName()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::phys::DiBoson<phys::Electron,phys::Electron>*)0x0)->GetImplFileName();
}

//______________________________________________________________________________
template <> int DiBoson<phys::Electron,phys::Electron>::ImplFileLine()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::phys::DiBoson<phys::Electron,phys::Electron>*)0x0)->GetImplFileLine();
}

//______________________________________________________________________________
template <> TClass *DiBoson<phys::Electron,phys::Electron>::Dictionary()
{
   fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::phys::DiBoson<phys::Electron,phys::Electron>*)0x0)->GetClass();
   return fgIsA;
}

//______________________________________________________________________________
template <> TClass *DiBoson<phys::Electron,phys::Electron>::Class()
{
   if (!fgIsA.load()) { R__LOCKGUARD(gInterpreterMutex); fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::phys::DiBoson<phys::Electron,phys::Electron>*)0x0)->GetClass(); }
   return fgIsA;
}

} // namespace phys
namespace phys {
//______________________________________________________________________________
template <> atomic_TClass_ptr DiBoson<phys::Electron,phys::Lepton>::fgIsA(0);  // static to hold class pointer

//______________________________________________________________________________
template <> const char *DiBoson<phys::Electron,phys::Lepton>::Class_Name()
{
   return "phys::DiBoson<phys::Electron,phys::Lepton>";
}

//______________________________________________________________________________
template <> const char *DiBoson<phys::Electron,phys::Lepton>::ImplFileName()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::phys::DiBoson<phys::Electron,phys::Lepton>*)0x0)->GetImplFileName();
}

//______________________________________________________________________________
template <> int DiBoson<phys::Electron,phys::Lepton>::ImplFileLine()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::phys::DiBoson<phys::Electron,phys::Lepton>*)0x0)->GetImplFileLine();
}

//______________________________________________________________________________
template <> TClass *DiBoson<phys::Electron,phys::Lepton>::Dictionary()
{
   fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::phys::DiBoson<phys::Electron,phys::Lepton>*)0x0)->GetClass();
   return fgIsA;
}

//______________________________________________________________________________
template <> TClass *DiBoson<phys::Electron,phys::Lepton>::Class()
{
   if (!fgIsA.load()) { R__LOCKGUARD(gInterpreterMutex); fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::phys::DiBoson<phys::Electron,phys::Lepton>*)0x0)->GetClass(); }
   return fgIsA;
}

} // namespace phys
namespace phys {
//______________________________________________________________________________
template <> atomic_TClass_ptr DiBoson<phys::Lepton,phys::Electron>::fgIsA(0);  // static to hold class pointer

//______________________________________________________________________________
template <> const char *DiBoson<phys::Lepton,phys::Electron>::Class_Name()
{
   return "phys::DiBoson<phys::Lepton,phys::Electron>";
}

//______________________________________________________________________________
template <> const char *DiBoson<phys::Lepton,phys::Electron>::ImplFileName()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::phys::DiBoson<phys::Lepton,phys::Electron>*)0x0)->GetImplFileName();
}

//______________________________________________________________________________
template <> int DiBoson<phys::Lepton,phys::Electron>::ImplFileLine()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::phys::DiBoson<phys::Lepton,phys::Electron>*)0x0)->GetImplFileLine();
}

//______________________________________________________________________________
template <> TClass *DiBoson<phys::Lepton,phys::Electron>::Dictionary()
{
   fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::phys::DiBoson<phys::Lepton,phys::Electron>*)0x0)->GetClass();
   return fgIsA;
}

//______________________________________________________________________________
template <> TClass *DiBoson<phys::Lepton,phys::Electron>::Class()
{
   if (!fgIsA.load()) { R__LOCKGUARD(gInterpreterMutex); fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::phys::DiBoson<phys::Lepton,phys::Electron>*)0x0)->GetClass(); }
   return fgIsA;
}

} // namespace phys
namespace phys {
//______________________________________________________________________________
template <> atomic_TClass_ptr DiBoson<phys::Particle,phys::Particle>::fgIsA(0);  // static to hold class pointer

//______________________________________________________________________________
template <> const char *DiBoson<phys::Particle,phys::Particle>::Class_Name()
{
   return "phys::DiBoson<phys::Particle,phys::Particle>";
}

//______________________________________________________________________________
template <> const char *DiBoson<phys::Particle,phys::Particle>::ImplFileName()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::phys::DiBoson<phys::Particle,phys::Particle>*)0x0)->GetImplFileName();
}

//______________________________________________________________________________
template <> int DiBoson<phys::Particle,phys::Particle>::ImplFileLine()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::phys::DiBoson<phys::Particle,phys::Particle>*)0x0)->GetImplFileLine();
}

//______________________________________________________________________________
template <> TClass *DiBoson<phys::Particle,phys::Particle>::Dictionary()
{
   fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::phys::DiBoson<phys::Particle,phys::Particle>*)0x0)->GetClass();
   return fgIsA;
}

//______________________________________________________________________________
template <> TClass *DiBoson<phys::Particle,phys::Particle>::Class()
{
   if (!fgIsA.load()) { R__LOCKGUARD(gInterpreterMutex); fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::phys::DiBoson<phys::Particle,phys::Particle>*)0x0)->GetClass(); }
   return fgIsA;
}

} // namespace phys
namespace phys {
//______________________________________________________________________________
atomic_TClass_ptr GenEventWeights::fgIsA(0);  // static to hold class pointer

//______________________________________________________________________________
const char *GenEventWeights::Class_Name()
{
   return "phys::GenEventWeights";
}

//______________________________________________________________________________
const char *GenEventWeights::ImplFileName()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::phys::GenEventWeights*)0x0)->GetImplFileName();
}

//______________________________________________________________________________
int GenEventWeights::ImplFileLine()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::phys::GenEventWeights*)0x0)->GetImplFileLine();
}

//______________________________________________________________________________
TClass *GenEventWeights::Dictionary()
{
   fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::phys::GenEventWeights*)0x0)->GetClass();
   return fgIsA;
}

//______________________________________________________________________________
TClass *GenEventWeights::Class()
{
   if (!fgIsA.load()) { R__LOCKGUARD(gInterpreterMutex); fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::phys::GenEventWeights*)0x0)->GetClass(); }
   return fgIsA;
}

} // namespace phys
namespace phys {
//______________________________________________________________________________
atomic_TClass_ptr MELA::fgIsA(0);  // static to hold class pointer

//______________________________________________________________________________
const char *MELA::Class_Name()
{
   return "phys::MELA";
}

//______________________________________________________________________________
const char *MELA::ImplFileName()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::phys::MELA*)0x0)->GetImplFileName();
}

//______________________________________________________________________________
int MELA::ImplFileLine()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::phys::MELA*)0x0)->GetImplFileLine();
}

//______________________________________________________________________________
TClass *MELA::Dictionary()
{
   fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::phys::MELA*)0x0)->GetClass();
   return fgIsA;
}

//______________________________________________________________________________
TClass *MELA::Class()
{
   if (!fgIsA.load()) { R__LOCKGUARD(gInterpreterMutex); fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::phys::MELA*)0x0)->GetClass(); }
   return fgIsA;
}

} // namespace phys
namespace ROOT {
   // Wrappers around operator new
   static void *new_pairlEphyscLcLBosonlEphyscLcLLeptongRcOphyscLcLLeptongR(void *p) {
      return  p ? ::new((::ROOT::Internal::TOperatorNewHelper*)p) pair<phys::Boson<phys::Lepton>,phys::Lepton> : new pair<phys::Boson<phys::Lepton>,phys::Lepton>;
   }
   static void *newArray_pairlEphyscLcLBosonlEphyscLcLLeptongRcOphyscLcLLeptongR(Long_t nElements, void *p) {
      return p ? ::new((::ROOT::Internal::TOperatorNewHelper*)p) pair<phys::Boson<phys::Lepton>,phys::Lepton>[nElements] : new pair<phys::Boson<phys::Lepton>,phys::Lepton>[nElements];
   }
   // Wrapper around operator delete
   static void delete_pairlEphyscLcLBosonlEphyscLcLLeptongRcOphyscLcLLeptongR(void *p) {
      delete ((pair<phys::Boson<phys::Lepton>,phys::Lepton>*)p);
   }
   static void deleteArray_pairlEphyscLcLBosonlEphyscLcLLeptongRcOphyscLcLLeptongR(void *p) {
      delete [] ((pair<phys::Boson<phys::Lepton>,phys::Lepton>*)p);
   }
   static void destruct_pairlEphyscLcLBosonlEphyscLcLLeptongRcOphyscLcLLeptongR(void *p) {
      typedef pair<phys::Boson<phys::Lepton>,phys::Lepton> current_t;
      ((current_t*)p)->~current_t();
   }
} // end of namespace ROOT for class pair<phys::Boson<phys::Lepton>,phys::Lepton>

namespace phys {
//______________________________________________________________________________
void Particle::Streamer(TBuffer &R__b)
{
   // Stream an object of class phys::Particle.

   if (R__b.IsReading()) {
      R__b.ReadClassBuffer(phys::Particle::Class(),this);
   } else {
      R__b.WriteClassBuffer(phys::Particle::Class(),this);
   }
}

} // namespace phys
namespace ROOT {
   // Wrappers around operator new
   static void *new_physcLcLParticle(void *p) {
      return  p ? new(p) ::phys::Particle : new ::phys::Particle;
   }
   static void *newArray_physcLcLParticle(Long_t nElements, void *p) {
      return p ? new(p) ::phys::Particle[nElements] : new ::phys::Particle[nElements];
   }
   // Wrapper around operator delete
   static void delete_physcLcLParticle(void *p) {
      delete ((::phys::Particle*)p);
   }
   static void deleteArray_physcLcLParticle(void *p) {
      delete [] ((::phys::Particle*)p);
   }
   static void destruct_physcLcLParticle(void *p) {
      typedef ::phys::Particle current_t;
      ((current_t*)p)->~current_t();
   }
} // end of namespace ROOT for class ::phys::Particle

namespace phys {
//______________________________________________________________________________
void Jet::Streamer(TBuffer &R__b)
{
   // Stream an object of class phys::Jet.

   if (R__b.IsReading()) {
      R__b.ReadClassBuffer(phys::Jet::Class(),this);
   } else {
      R__b.WriteClassBuffer(phys::Jet::Class(),this);
   }
}

} // namespace phys
namespace ROOT {
   // Wrappers around operator new
   static void *new_physcLcLJet(void *p) {
      return  p ? new(p) ::phys::Jet : new ::phys::Jet;
   }
   static void *newArray_physcLcLJet(Long_t nElements, void *p) {
      return p ? new(p) ::phys::Jet[nElements] : new ::phys::Jet[nElements];
   }
   // Wrapper around operator delete
   static void delete_physcLcLJet(void *p) {
      delete ((::phys::Jet*)p);
   }
   static void deleteArray_physcLcLJet(void *p) {
      delete [] ((::phys::Jet*)p);
   }
   static void destruct_physcLcLJet(void *p) {
      typedef ::phys::Jet current_t;
      ((current_t*)p)->~current_t();
   }
} // end of namespace ROOT for class ::phys::Jet

namespace phys {
//______________________________________________________________________________
void Lepton::Streamer(TBuffer &R__b)
{
   // Stream an object of class phys::Lepton.

   if (R__b.IsReading()) {
      R__b.ReadClassBuffer(phys::Lepton::Class(),this);
   } else {
      R__b.WriteClassBuffer(phys::Lepton::Class(),this);
   }
}

} // namespace phys
namespace ROOT {
   // Wrappers around operator new
   static void *new_physcLcLLepton(void *p) {
      return  p ? new(p) ::phys::Lepton : new ::phys::Lepton;
   }
   static void *newArray_physcLcLLepton(Long_t nElements, void *p) {
      return p ? new(p) ::phys::Lepton[nElements] : new ::phys::Lepton[nElements];
   }
   // Wrapper around operator delete
   static void delete_physcLcLLepton(void *p) {
      delete ((::phys::Lepton*)p);
   }
   static void deleteArray_physcLcLLepton(void *p) {
      delete [] ((::phys::Lepton*)p);
   }
   static void destruct_physcLcLLepton(void *p) {
      typedef ::phys::Lepton current_t;
      ((current_t*)p)->~current_t();
   }
} // end of namespace ROOT for class ::phys::Lepton

namespace phys {
//______________________________________________________________________________
void Electron::Streamer(TBuffer &R__b)
{
   // Stream an object of class phys::Electron.

   if (R__b.IsReading()) {
      R__b.ReadClassBuffer(phys::Electron::Class(),this);
   } else {
      R__b.WriteClassBuffer(phys::Electron::Class(),this);
   }
}

} // namespace phys
namespace ROOT {
   // Wrappers around operator new
   static void *new_physcLcLElectron(void *p) {
      return  p ? new(p) ::phys::Electron : new ::phys::Electron;
   }
   static void *newArray_physcLcLElectron(Long_t nElements, void *p) {
      return p ? new(p) ::phys::Electron[nElements] : new ::phys::Electron[nElements];
   }
   // Wrapper around operator delete
   static void delete_physcLcLElectron(void *p) {
      delete ((::phys::Electron*)p);
   }
   static void deleteArray_physcLcLElectron(void *p) {
      delete [] ((::phys::Electron*)p);
   }
   static void destruct_physcLcLElectron(void *p) {
      typedef ::phys::Electron current_t;
      ((current_t*)p)->~current_t();
   }
} // end of namespace ROOT for class ::phys::Electron

namespace phys {
//______________________________________________________________________________
template <> void Boson<phys::Particle>::Streamer(TBuffer &R__b)
{
   // Stream an object of class phys::Boson<phys::Particle>.

   if (R__b.IsReading()) {
      R__b.ReadClassBuffer(phys::Boson<phys::Particle>::Class(),this);
   } else {
      R__b.WriteClassBuffer(phys::Boson<phys::Particle>::Class(),this);
   }
}

} // namespace phys
namespace ROOT {
   // Wrappers around operator new
   static void *new_physcLcLBosonlEphyscLcLParticlegR(void *p) {
      return  p ? new(p) ::phys::Boson<phys::Particle> : new ::phys::Boson<phys::Particle>;
   }
   static void *newArray_physcLcLBosonlEphyscLcLParticlegR(Long_t nElements, void *p) {
      return p ? new(p) ::phys::Boson<phys::Particle>[nElements] : new ::phys::Boson<phys::Particle>[nElements];
   }
   // Wrapper around operator delete
   static void delete_physcLcLBosonlEphyscLcLParticlegR(void *p) {
      delete ((::phys::Boson<phys::Particle>*)p);
   }
   static void deleteArray_physcLcLBosonlEphyscLcLParticlegR(void *p) {
      delete [] ((::phys::Boson<phys::Particle>*)p);
   }
   static void destruct_physcLcLBosonlEphyscLcLParticlegR(void *p) {
      typedef ::phys::Boson<phys::Particle> current_t;
      ((current_t*)p)->~current_t();
   }
} // end of namespace ROOT for class ::phys::Boson<phys::Particle>

namespace phys {
//______________________________________________________________________________
template <> void Boson<phys::Lepton>::Streamer(TBuffer &R__b)
{
   // Stream an object of class phys::Boson<phys::Lepton>.

   if (R__b.IsReading()) {
      R__b.ReadClassBuffer(phys::Boson<phys::Lepton>::Class(),this);
   } else {
      R__b.WriteClassBuffer(phys::Boson<phys::Lepton>::Class(),this);
   }
}

} // namespace phys
namespace ROOT {
   // Wrappers around operator new
   static void *new_physcLcLBosonlEphyscLcLLeptongR(void *p) {
      return  p ? new(p) ::phys::Boson<phys::Lepton> : new ::phys::Boson<phys::Lepton>;
   }
   static void *newArray_physcLcLBosonlEphyscLcLLeptongR(Long_t nElements, void *p) {
      return p ? new(p) ::phys::Boson<phys::Lepton>[nElements] : new ::phys::Boson<phys::Lepton>[nElements];
   }
   // Wrapper around operator delete
   static void delete_physcLcLBosonlEphyscLcLLeptongR(void *p) {
      delete ((::phys::Boson<phys::Lepton>*)p);
   }
   static void deleteArray_physcLcLBosonlEphyscLcLLeptongR(void *p) {
      delete [] ((::phys::Boson<phys::Lepton>*)p);
   }
   static void destruct_physcLcLBosonlEphyscLcLLeptongR(void *p) {
      typedef ::phys::Boson<phys::Lepton> current_t;
      ((current_t*)p)->~current_t();
   }
} // end of namespace ROOT for class ::phys::Boson<phys::Lepton>

namespace phys {
//______________________________________________________________________________
template <> void Boson<phys::Jet>::Streamer(TBuffer &R__b)
{
   // Stream an object of class phys::Boson<phys::Jet>.

   if (R__b.IsReading()) {
      R__b.ReadClassBuffer(phys::Boson<phys::Jet>::Class(),this);
   } else {
      R__b.WriteClassBuffer(phys::Boson<phys::Jet>::Class(),this);
   }
}

} // namespace phys
namespace ROOT {
   // Wrappers around operator new
   static void *new_physcLcLBosonlEphyscLcLJetgR(void *p) {
      return  p ? new(p) ::phys::Boson<phys::Jet> : new ::phys::Boson<phys::Jet>;
   }
   static void *newArray_physcLcLBosonlEphyscLcLJetgR(Long_t nElements, void *p) {
      return p ? new(p) ::phys::Boson<phys::Jet>[nElements] : new ::phys::Boson<phys::Jet>[nElements];
   }
   // Wrapper around operator delete
   static void delete_physcLcLBosonlEphyscLcLJetgR(void *p) {
      delete ((::phys::Boson<phys::Jet>*)p);
   }
   static void deleteArray_physcLcLBosonlEphyscLcLJetgR(void *p) {
      delete [] ((::phys::Boson<phys::Jet>*)p);
   }
   static void destruct_physcLcLBosonlEphyscLcLJetgR(void *p) {
      typedef ::phys::Boson<phys::Jet> current_t;
      ((current_t*)p)->~current_t();
   }
} // end of namespace ROOT for class ::phys::Boson<phys::Jet>

namespace phys {
//______________________________________________________________________________
template <> void Boson<phys::Electron>::Streamer(TBuffer &R__b)
{
   // Stream an object of class phys::Boson<phys::Electron>.

   if (R__b.IsReading()) {
      R__b.ReadClassBuffer(phys::Boson<phys::Electron>::Class(),this);
   } else {
      R__b.WriteClassBuffer(phys::Boson<phys::Electron>::Class(),this);
   }
}

} // namespace phys
namespace ROOT {
   // Wrappers around operator new
   static void *new_physcLcLBosonlEphyscLcLElectrongR(void *p) {
      return  p ? new(p) ::phys::Boson<phys::Electron> : new ::phys::Boson<phys::Electron>;
   }
   static void *newArray_physcLcLBosonlEphyscLcLElectrongR(Long_t nElements, void *p) {
      return p ? new(p) ::phys::Boson<phys::Electron>[nElements] : new ::phys::Boson<phys::Electron>[nElements];
   }
   // Wrapper around operator delete
   static void delete_physcLcLBosonlEphyscLcLElectrongR(void *p) {
      delete ((::phys::Boson<phys::Electron>*)p);
   }
   static void deleteArray_physcLcLBosonlEphyscLcLElectrongR(void *p) {
      delete [] ((::phys::Boson<phys::Electron>*)p);
   }
   static void destruct_physcLcLBosonlEphyscLcLElectrongR(void *p) {
      typedef ::phys::Boson<phys::Electron> current_t;
      ((current_t*)p)->~current_t();
   }
} // end of namespace ROOT for class ::phys::Boson<phys::Electron>

namespace phys {
//______________________________________________________________________________
template <> void DiBoson<phys::Lepton,phys::Lepton>::Streamer(TBuffer &R__b)
{
   // Stream an object of class phys::DiBoson<phys::Lepton,phys::Lepton>.

   if (R__b.IsReading()) {
      R__b.ReadClassBuffer(phys::DiBoson<phys::Lepton,phys::Lepton>::Class(),this);
   } else {
      R__b.WriteClassBuffer(phys::DiBoson<phys::Lepton,phys::Lepton>::Class(),this);
   }
}

} // namespace phys
namespace ROOT {
   // Wrappers around operator new
   static void *new_physcLcLDiBosonlEphyscLcLLeptoncOphyscLcLLeptongR(void *p) {
      return  p ? new(p) ::phys::DiBoson<phys::Lepton,phys::Lepton> : new ::phys::DiBoson<phys::Lepton,phys::Lepton>;
   }
   static void *newArray_physcLcLDiBosonlEphyscLcLLeptoncOphyscLcLLeptongR(Long_t nElements, void *p) {
      return p ? new(p) ::phys::DiBoson<phys::Lepton,phys::Lepton>[nElements] : new ::phys::DiBoson<phys::Lepton,phys::Lepton>[nElements];
   }
   // Wrapper around operator delete
   static void delete_physcLcLDiBosonlEphyscLcLLeptoncOphyscLcLLeptongR(void *p) {
      delete ((::phys::DiBoson<phys::Lepton,phys::Lepton>*)p);
   }
   static void deleteArray_physcLcLDiBosonlEphyscLcLLeptoncOphyscLcLLeptongR(void *p) {
      delete [] ((::phys::DiBoson<phys::Lepton,phys::Lepton>*)p);
   }
   static void destruct_physcLcLDiBosonlEphyscLcLLeptoncOphyscLcLLeptongR(void *p) {
      typedef ::phys::DiBoson<phys::Lepton,phys::Lepton> current_t;
      ((current_t*)p)->~current_t();
   }
} // end of namespace ROOT for class ::phys::DiBoson<phys::Lepton,phys::Lepton>

namespace phys {
//______________________________________________________________________________
template <> void DiBoson<phys::Electron,phys::Electron>::Streamer(TBuffer &R__b)
{
   // Stream an object of class phys::DiBoson<phys::Electron,phys::Electron>.

   if (R__b.IsReading()) {
      R__b.ReadClassBuffer(phys::DiBoson<phys::Electron,phys::Electron>::Class(),this);
   } else {
      R__b.WriteClassBuffer(phys::DiBoson<phys::Electron,phys::Electron>::Class(),this);
   }
}

} // namespace phys
namespace ROOT {
   // Wrappers around operator new
   static void *new_physcLcLDiBosonlEphyscLcLElectroncOphyscLcLElectrongR(void *p) {
      return  p ? new(p) ::phys::DiBoson<phys::Electron,phys::Electron> : new ::phys::DiBoson<phys::Electron,phys::Electron>;
   }
   static void *newArray_physcLcLDiBosonlEphyscLcLElectroncOphyscLcLElectrongR(Long_t nElements, void *p) {
      return p ? new(p) ::phys::DiBoson<phys::Electron,phys::Electron>[nElements] : new ::phys::DiBoson<phys::Electron,phys::Electron>[nElements];
   }
   // Wrapper around operator delete
   static void delete_physcLcLDiBosonlEphyscLcLElectroncOphyscLcLElectrongR(void *p) {
      delete ((::phys::DiBoson<phys::Electron,phys::Electron>*)p);
   }
   static void deleteArray_physcLcLDiBosonlEphyscLcLElectroncOphyscLcLElectrongR(void *p) {
      delete [] ((::phys::DiBoson<phys::Electron,phys::Electron>*)p);
   }
   static void destruct_physcLcLDiBosonlEphyscLcLElectroncOphyscLcLElectrongR(void *p) {
      typedef ::phys::DiBoson<phys::Electron,phys::Electron> current_t;
      ((current_t*)p)->~current_t();
   }
} // end of namespace ROOT for class ::phys::DiBoson<phys::Electron,phys::Electron>

namespace phys {
//______________________________________________________________________________
template <> void DiBoson<phys::Electron,phys::Lepton>::Streamer(TBuffer &R__b)
{
   // Stream an object of class phys::DiBoson<phys::Electron,phys::Lepton>.

   if (R__b.IsReading()) {
      R__b.ReadClassBuffer(phys::DiBoson<phys::Electron,phys::Lepton>::Class(),this);
   } else {
      R__b.WriteClassBuffer(phys::DiBoson<phys::Electron,phys::Lepton>::Class(),this);
   }
}

} // namespace phys
namespace ROOT {
   // Wrappers around operator new
   static void *new_physcLcLDiBosonlEphyscLcLElectroncOphyscLcLLeptongR(void *p) {
      return  p ? new(p) ::phys::DiBoson<phys::Electron,phys::Lepton> : new ::phys::DiBoson<phys::Electron,phys::Lepton>;
   }
   static void *newArray_physcLcLDiBosonlEphyscLcLElectroncOphyscLcLLeptongR(Long_t nElements, void *p) {
      return p ? new(p) ::phys::DiBoson<phys::Electron,phys::Lepton>[nElements] : new ::phys::DiBoson<phys::Electron,phys::Lepton>[nElements];
   }
   // Wrapper around operator delete
   static void delete_physcLcLDiBosonlEphyscLcLElectroncOphyscLcLLeptongR(void *p) {
      delete ((::phys::DiBoson<phys::Electron,phys::Lepton>*)p);
   }
   static void deleteArray_physcLcLDiBosonlEphyscLcLElectroncOphyscLcLLeptongR(void *p) {
      delete [] ((::phys::DiBoson<phys::Electron,phys::Lepton>*)p);
   }
   static void destruct_physcLcLDiBosonlEphyscLcLElectroncOphyscLcLLeptongR(void *p) {
      typedef ::phys::DiBoson<phys::Electron,phys::Lepton> current_t;
      ((current_t*)p)->~current_t();
   }
} // end of namespace ROOT for class ::phys::DiBoson<phys::Electron,phys::Lepton>

namespace phys {
//______________________________________________________________________________
template <> void DiBoson<phys::Lepton,phys::Electron>::Streamer(TBuffer &R__b)
{
   // Stream an object of class phys::DiBoson<phys::Lepton,phys::Electron>.

   if (R__b.IsReading()) {
      R__b.ReadClassBuffer(phys::DiBoson<phys::Lepton,phys::Electron>::Class(),this);
   } else {
      R__b.WriteClassBuffer(phys::DiBoson<phys::Lepton,phys::Electron>::Class(),this);
   }
}

} // namespace phys
namespace ROOT {
   // Wrappers around operator new
   static void *new_physcLcLDiBosonlEphyscLcLLeptoncOphyscLcLElectrongR(void *p) {
      return  p ? new(p) ::phys::DiBoson<phys::Lepton,phys::Electron> : new ::phys::DiBoson<phys::Lepton,phys::Electron>;
   }
   static void *newArray_physcLcLDiBosonlEphyscLcLLeptoncOphyscLcLElectrongR(Long_t nElements, void *p) {
      return p ? new(p) ::phys::DiBoson<phys::Lepton,phys::Electron>[nElements] : new ::phys::DiBoson<phys::Lepton,phys::Electron>[nElements];
   }
   // Wrapper around operator delete
   static void delete_physcLcLDiBosonlEphyscLcLLeptoncOphyscLcLElectrongR(void *p) {
      delete ((::phys::DiBoson<phys::Lepton,phys::Electron>*)p);
   }
   static void deleteArray_physcLcLDiBosonlEphyscLcLLeptoncOphyscLcLElectrongR(void *p) {
      delete [] ((::phys::DiBoson<phys::Lepton,phys::Electron>*)p);
   }
   static void destruct_physcLcLDiBosonlEphyscLcLLeptoncOphyscLcLElectrongR(void *p) {
      typedef ::phys::DiBoson<phys::Lepton,phys::Electron> current_t;
      ((current_t*)p)->~current_t();
   }
} // end of namespace ROOT for class ::phys::DiBoson<phys::Lepton,phys::Electron>

namespace phys {
//______________________________________________________________________________
template <> void DiBoson<phys::Particle,phys::Particle>::Streamer(TBuffer &R__b)
{
   // Stream an object of class phys::DiBoson<phys::Particle,phys::Particle>.

   if (R__b.IsReading()) {
      R__b.ReadClassBuffer(phys::DiBoson<phys::Particle,phys::Particle>::Class(),this);
   } else {
      R__b.WriteClassBuffer(phys::DiBoson<phys::Particle,phys::Particle>::Class(),this);
   }
}

} // namespace phys
namespace ROOT {
   // Wrappers around operator new
   static void *new_physcLcLDiBosonlEphyscLcLParticlecOphyscLcLParticlegR(void *p) {
      return  p ? new(p) ::phys::DiBoson<phys::Particle,phys::Particle> : new ::phys::DiBoson<phys::Particle,phys::Particle>;
   }
   static void *newArray_physcLcLDiBosonlEphyscLcLParticlecOphyscLcLParticlegR(Long_t nElements, void *p) {
      return p ? new(p) ::phys::DiBoson<phys::Particle,phys::Particle>[nElements] : new ::phys::DiBoson<phys::Particle,phys::Particle>[nElements];
   }
   // Wrapper around operator delete
   static void delete_physcLcLDiBosonlEphyscLcLParticlecOphyscLcLParticlegR(void *p) {
      delete ((::phys::DiBoson<phys::Particle,phys::Particle>*)p);
   }
   static void deleteArray_physcLcLDiBosonlEphyscLcLParticlecOphyscLcLParticlegR(void *p) {
      delete [] ((::phys::DiBoson<phys::Particle,phys::Particle>*)p);
   }
   static void destruct_physcLcLDiBosonlEphyscLcLParticlecOphyscLcLParticlegR(void *p) {
      typedef ::phys::DiBoson<phys::Particle,phys::Particle> current_t;
      ((current_t*)p)->~current_t();
   }
} // end of namespace ROOT for class ::phys::DiBoson<phys::Particle,phys::Particle>

namespace phys {
//______________________________________________________________________________
void GenEventWeights::Streamer(TBuffer &R__b)
{
   // Stream an object of class phys::GenEventWeights.

   if (R__b.IsReading()) {
      R__b.ReadClassBuffer(phys::GenEventWeights::Class(),this);
   } else {
      R__b.WriteClassBuffer(phys::GenEventWeights::Class(),this);
   }
}

} // namespace phys
namespace ROOT {
   // Wrappers around operator new
   static void *new_physcLcLGenEventWeights(void *p) {
      return  p ? new(p) ::phys::GenEventWeights : new ::phys::GenEventWeights;
   }
   static void *newArray_physcLcLGenEventWeights(Long_t nElements, void *p) {
      return p ? new(p) ::phys::GenEventWeights[nElements] : new ::phys::GenEventWeights[nElements];
   }
   // Wrapper around operator delete
   static void delete_physcLcLGenEventWeights(void *p) {
      delete ((::phys::GenEventWeights*)p);
   }
   static void deleteArray_physcLcLGenEventWeights(void *p) {
      delete [] ((::phys::GenEventWeights*)p);
   }
   static void destruct_physcLcLGenEventWeights(void *p) {
      typedef ::phys::GenEventWeights current_t;
      ((current_t*)p)->~current_t();
   }
} // end of namespace ROOT for class ::phys::GenEventWeights

namespace phys {
//______________________________________________________________________________
void MELA::Streamer(TBuffer &R__b)
{
   // Stream an object of class phys::MELA.

   if (R__b.IsReading()) {
      R__b.ReadClassBuffer(phys::MELA::Class(),this);
   } else {
      R__b.WriteClassBuffer(phys::MELA::Class(),this);
   }
}

} // namespace phys
namespace ROOT {
   // Wrappers around operator new
   static void *new_physcLcLMELA(void *p) {
      return  p ? new(p) ::phys::MELA : new ::phys::MELA;
   }
   static void *newArray_physcLcLMELA(Long_t nElements, void *p) {
      return p ? new(p) ::phys::MELA[nElements] : new ::phys::MELA[nElements];
   }
   // Wrapper around operator delete
   static void delete_physcLcLMELA(void *p) {
      delete ((::phys::MELA*)p);
   }
   static void deleteArray_physcLcLMELA(void *p) {
      delete [] ((::phys::MELA*)p);
   }
   static void destruct_physcLcLMELA(void *p) {
      typedef ::phys::MELA current_t;
      ((current_t*)p)->~current_t();
   }
} // end of namespace ROOT for class ::phys::MELA

namespace ROOT {
   static TClass *vectorlEphyscLcLParticlegR_Dictionary();
   static void vectorlEphyscLcLParticlegR_TClassManip(TClass*);
   static void *new_vectorlEphyscLcLParticlegR(void *p = 0);
   static void *newArray_vectorlEphyscLcLParticlegR(Long_t size, void *p);
   static void delete_vectorlEphyscLcLParticlegR(void *p);
   static void deleteArray_vectorlEphyscLcLParticlegR(void *p);
   static void destruct_vectorlEphyscLcLParticlegR(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const vector<phys::Particle>*)
   {
      vector<phys::Particle> *ptr = 0;
      static ::TVirtualIsAProxy* isa_proxy = new ::TIsAProxy(typeid(vector<phys::Particle>));
      static ::ROOT::TGenericClassInfo 
         instance("vector<phys::Particle>", -2, "vector", 214,
                  typeid(vector<phys::Particle>), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &vectorlEphyscLcLParticlegR_Dictionary, isa_proxy, 4,
                  sizeof(vector<phys::Particle>) );
      instance.SetNew(&new_vectorlEphyscLcLParticlegR);
      instance.SetNewArray(&newArray_vectorlEphyscLcLParticlegR);
      instance.SetDelete(&delete_vectorlEphyscLcLParticlegR);
      instance.SetDeleteArray(&deleteArray_vectorlEphyscLcLParticlegR);
      instance.SetDestructor(&destruct_vectorlEphyscLcLParticlegR);
      instance.AdoptCollectionProxyInfo(TCollectionProxyInfo::Generate(TCollectionProxyInfo::Pushback< vector<phys::Particle> >()));
      return &instance;
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal((const vector<phys::Particle>*)0x0); R__UseDummy(_R__UNIQUE_DICT_(Init));

   // Dictionary for non-ClassDef classes
   static TClass *vectorlEphyscLcLParticlegR_Dictionary() {
      TClass* theClass =::ROOT::GenerateInitInstanceLocal((const vector<phys::Particle>*)0x0)->GetClass();
      vectorlEphyscLcLParticlegR_TClassManip(theClass);
   return theClass;
   }

   static void vectorlEphyscLcLParticlegR_TClassManip(TClass* ){
   }

} // end of namespace ROOT

namespace ROOT {
   // Wrappers around operator new
   static void *new_vectorlEphyscLcLParticlegR(void *p) {
      return  p ? ::new((::ROOT::Internal::TOperatorNewHelper*)p) vector<phys::Particle> : new vector<phys::Particle>;
   }
   static void *newArray_vectorlEphyscLcLParticlegR(Long_t nElements, void *p) {
      return p ? ::new((::ROOT::Internal::TOperatorNewHelper*)p) vector<phys::Particle>[nElements] : new vector<phys::Particle>[nElements];
   }
   // Wrapper around operator delete
   static void delete_vectorlEphyscLcLParticlegR(void *p) {
      delete ((vector<phys::Particle>*)p);
   }
   static void deleteArray_vectorlEphyscLcLParticlegR(void *p) {
      delete [] ((vector<phys::Particle>*)p);
   }
   static void destruct_vectorlEphyscLcLParticlegR(void *p) {
      typedef vector<phys::Particle> current_t;
      ((current_t*)p)->~current_t();
   }
} // end of namespace ROOT for class vector<phys::Particle>

namespace ROOT {
   static TClass *vectorlEphyscLcLLeptongR_Dictionary();
   static void vectorlEphyscLcLLeptongR_TClassManip(TClass*);
   static void *new_vectorlEphyscLcLLeptongR(void *p = 0);
   static void *newArray_vectorlEphyscLcLLeptongR(Long_t size, void *p);
   static void delete_vectorlEphyscLcLLeptongR(void *p);
   static void deleteArray_vectorlEphyscLcLLeptongR(void *p);
   static void destruct_vectorlEphyscLcLLeptongR(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const vector<phys::Lepton>*)
   {
      vector<phys::Lepton> *ptr = 0;
      static ::TVirtualIsAProxy* isa_proxy = new ::TIsAProxy(typeid(vector<phys::Lepton>));
      static ::ROOT::TGenericClassInfo 
         instance("vector<phys::Lepton>", -2, "vector", 214,
                  typeid(vector<phys::Lepton>), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &vectorlEphyscLcLLeptongR_Dictionary, isa_proxy, 4,
                  sizeof(vector<phys::Lepton>) );
      instance.SetNew(&new_vectorlEphyscLcLLeptongR);
      instance.SetNewArray(&newArray_vectorlEphyscLcLLeptongR);
      instance.SetDelete(&delete_vectorlEphyscLcLLeptongR);
      instance.SetDeleteArray(&deleteArray_vectorlEphyscLcLLeptongR);
      instance.SetDestructor(&destruct_vectorlEphyscLcLLeptongR);
      instance.AdoptCollectionProxyInfo(TCollectionProxyInfo::Generate(TCollectionProxyInfo::Pushback< vector<phys::Lepton> >()));
      return &instance;
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal((const vector<phys::Lepton>*)0x0); R__UseDummy(_R__UNIQUE_DICT_(Init));

   // Dictionary for non-ClassDef classes
   static TClass *vectorlEphyscLcLLeptongR_Dictionary() {
      TClass* theClass =::ROOT::GenerateInitInstanceLocal((const vector<phys::Lepton>*)0x0)->GetClass();
      vectorlEphyscLcLLeptongR_TClassManip(theClass);
   return theClass;
   }

   static void vectorlEphyscLcLLeptongR_TClassManip(TClass* ){
   }

} // end of namespace ROOT

namespace ROOT {
   // Wrappers around operator new
   static void *new_vectorlEphyscLcLLeptongR(void *p) {
      return  p ? ::new((::ROOT::Internal::TOperatorNewHelper*)p) vector<phys::Lepton> : new vector<phys::Lepton>;
   }
   static void *newArray_vectorlEphyscLcLLeptongR(Long_t nElements, void *p) {
      return p ? ::new((::ROOT::Internal::TOperatorNewHelper*)p) vector<phys::Lepton>[nElements] : new vector<phys::Lepton>[nElements];
   }
   // Wrapper around operator delete
   static void delete_vectorlEphyscLcLLeptongR(void *p) {
      delete ((vector<phys::Lepton>*)p);
   }
   static void deleteArray_vectorlEphyscLcLLeptongR(void *p) {
      delete [] ((vector<phys::Lepton>*)p);
   }
   static void destruct_vectorlEphyscLcLLeptongR(void *p) {
      typedef vector<phys::Lepton> current_t;
      ((current_t*)p)->~current_t();
   }
} // end of namespace ROOT for class vector<phys::Lepton>

namespace ROOT {
   static TClass *vectorlEphyscLcLJetgR_Dictionary();
   static void vectorlEphyscLcLJetgR_TClassManip(TClass*);
   static void *new_vectorlEphyscLcLJetgR(void *p = 0);
   static void *newArray_vectorlEphyscLcLJetgR(Long_t size, void *p);
   static void delete_vectorlEphyscLcLJetgR(void *p);
   static void deleteArray_vectorlEphyscLcLJetgR(void *p);
   static void destruct_vectorlEphyscLcLJetgR(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const vector<phys::Jet>*)
   {
      vector<phys::Jet> *ptr = 0;
      static ::TVirtualIsAProxy* isa_proxy = new ::TIsAProxy(typeid(vector<phys::Jet>));
      static ::ROOT::TGenericClassInfo 
         instance("vector<phys::Jet>", -2, "vector", 214,
                  typeid(vector<phys::Jet>), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &vectorlEphyscLcLJetgR_Dictionary, isa_proxy, 4,
                  sizeof(vector<phys::Jet>) );
      instance.SetNew(&new_vectorlEphyscLcLJetgR);
      instance.SetNewArray(&newArray_vectorlEphyscLcLJetgR);
      instance.SetDelete(&delete_vectorlEphyscLcLJetgR);
      instance.SetDeleteArray(&deleteArray_vectorlEphyscLcLJetgR);
      instance.SetDestructor(&destruct_vectorlEphyscLcLJetgR);
      instance.AdoptCollectionProxyInfo(TCollectionProxyInfo::Generate(TCollectionProxyInfo::Pushback< vector<phys::Jet> >()));
      return &instance;
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal((const vector<phys::Jet>*)0x0); R__UseDummy(_R__UNIQUE_DICT_(Init));

   // Dictionary for non-ClassDef classes
   static TClass *vectorlEphyscLcLJetgR_Dictionary() {
      TClass* theClass =::ROOT::GenerateInitInstanceLocal((const vector<phys::Jet>*)0x0)->GetClass();
      vectorlEphyscLcLJetgR_TClassManip(theClass);
   return theClass;
   }

   static void vectorlEphyscLcLJetgR_TClassManip(TClass* ){
   }

} // end of namespace ROOT

namespace ROOT {
   // Wrappers around operator new
   static void *new_vectorlEphyscLcLJetgR(void *p) {
      return  p ? ::new((::ROOT::Internal::TOperatorNewHelper*)p) vector<phys::Jet> : new vector<phys::Jet>;
   }
   static void *newArray_vectorlEphyscLcLJetgR(Long_t nElements, void *p) {
      return p ? ::new((::ROOT::Internal::TOperatorNewHelper*)p) vector<phys::Jet>[nElements] : new vector<phys::Jet>[nElements];
   }
   // Wrapper around operator delete
   static void delete_vectorlEphyscLcLJetgR(void *p) {
      delete ((vector<phys::Jet>*)p);
   }
   static void deleteArray_vectorlEphyscLcLJetgR(void *p) {
      delete [] ((vector<phys::Jet>*)p);
   }
   static void destruct_vectorlEphyscLcLJetgR(void *p) {
      typedef vector<phys::Jet> current_t;
      ((current_t*)p)->~current_t();
   }
} // end of namespace ROOT for class vector<phys::Jet>

namespace ROOT {
   static TClass *vectorlEphyscLcLElectrongR_Dictionary();
   static void vectorlEphyscLcLElectrongR_TClassManip(TClass*);
   static void *new_vectorlEphyscLcLElectrongR(void *p = 0);
   static void *newArray_vectorlEphyscLcLElectrongR(Long_t size, void *p);
   static void delete_vectorlEphyscLcLElectrongR(void *p);
   static void deleteArray_vectorlEphyscLcLElectrongR(void *p);
   static void destruct_vectorlEphyscLcLElectrongR(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const vector<phys::Electron>*)
   {
      vector<phys::Electron> *ptr = 0;
      static ::TVirtualIsAProxy* isa_proxy = new ::TIsAProxy(typeid(vector<phys::Electron>));
      static ::ROOT::TGenericClassInfo 
         instance("vector<phys::Electron>", -2, "vector", 214,
                  typeid(vector<phys::Electron>), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &vectorlEphyscLcLElectrongR_Dictionary, isa_proxy, 4,
                  sizeof(vector<phys::Electron>) );
      instance.SetNew(&new_vectorlEphyscLcLElectrongR);
      instance.SetNewArray(&newArray_vectorlEphyscLcLElectrongR);
      instance.SetDelete(&delete_vectorlEphyscLcLElectrongR);
      instance.SetDeleteArray(&deleteArray_vectorlEphyscLcLElectrongR);
      instance.SetDestructor(&destruct_vectorlEphyscLcLElectrongR);
      instance.AdoptCollectionProxyInfo(TCollectionProxyInfo::Generate(TCollectionProxyInfo::Pushback< vector<phys::Electron> >()));
      return &instance;
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal((const vector<phys::Electron>*)0x0); R__UseDummy(_R__UNIQUE_DICT_(Init));

   // Dictionary for non-ClassDef classes
   static TClass *vectorlEphyscLcLElectrongR_Dictionary() {
      TClass* theClass =::ROOT::GenerateInitInstanceLocal((const vector<phys::Electron>*)0x0)->GetClass();
      vectorlEphyscLcLElectrongR_TClassManip(theClass);
   return theClass;
   }

   static void vectorlEphyscLcLElectrongR_TClassManip(TClass* ){
   }

} // end of namespace ROOT

namespace ROOT {
   // Wrappers around operator new
   static void *new_vectorlEphyscLcLElectrongR(void *p) {
      return  p ? ::new((::ROOT::Internal::TOperatorNewHelper*)p) vector<phys::Electron> : new vector<phys::Electron>;
   }
   static void *newArray_vectorlEphyscLcLElectrongR(Long_t nElements, void *p) {
      return p ? ::new((::ROOT::Internal::TOperatorNewHelper*)p) vector<phys::Electron>[nElements] : new vector<phys::Electron>[nElements];
   }
   // Wrapper around operator delete
   static void delete_vectorlEphyscLcLElectrongR(void *p) {
      delete ((vector<phys::Electron>*)p);
   }
   static void deleteArray_vectorlEphyscLcLElectrongR(void *p) {
      delete [] ((vector<phys::Electron>*)p);
   }
   static void destruct_vectorlEphyscLcLElectrongR(void *p) {
      typedef vector<phys::Electron> current_t;
      ((current_t*)p)->~current_t();
   }
} // end of namespace ROOT for class vector<phys::Electron>

namespace ROOT {
   static TClass *vectorlEphyscLcLDiBosonlEphyscLcLParticlecOphyscLcLParticlegRsPgR_Dictionary();
   static void vectorlEphyscLcLDiBosonlEphyscLcLParticlecOphyscLcLParticlegRsPgR_TClassManip(TClass*);
   static void *new_vectorlEphyscLcLDiBosonlEphyscLcLParticlecOphyscLcLParticlegRsPgR(void *p = 0);
   static void *newArray_vectorlEphyscLcLDiBosonlEphyscLcLParticlecOphyscLcLParticlegRsPgR(Long_t size, void *p);
   static void delete_vectorlEphyscLcLDiBosonlEphyscLcLParticlecOphyscLcLParticlegRsPgR(void *p);
   static void deleteArray_vectorlEphyscLcLDiBosonlEphyscLcLParticlecOphyscLcLParticlegRsPgR(void *p);
   static void destruct_vectorlEphyscLcLDiBosonlEphyscLcLParticlecOphyscLcLParticlegRsPgR(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const vector<phys::DiBoson<phys::Particle,phys::Particle> >*)
   {
      vector<phys::DiBoson<phys::Particle,phys::Particle> > *ptr = 0;
      static ::TVirtualIsAProxy* isa_proxy = new ::TIsAProxy(typeid(vector<phys::DiBoson<phys::Particle,phys::Particle> >));
      static ::ROOT::TGenericClassInfo 
         instance("vector<phys::DiBoson<phys::Particle,phys::Particle> >", -2, "vector", 214,
                  typeid(vector<phys::DiBoson<phys::Particle,phys::Particle> >), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &vectorlEphyscLcLDiBosonlEphyscLcLParticlecOphyscLcLParticlegRsPgR_Dictionary, isa_proxy, 4,
                  sizeof(vector<phys::DiBoson<phys::Particle,phys::Particle> >) );
      instance.SetNew(&new_vectorlEphyscLcLDiBosonlEphyscLcLParticlecOphyscLcLParticlegRsPgR);
      instance.SetNewArray(&newArray_vectorlEphyscLcLDiBosonlEphyscLcLParticlecOphyscLcLParticlegRsPgR);
      instance.SetDelete(&delete_vectorlEphyscLcLDiBosonlEphyscLcLParticlecOphyscLcLParticlegRsPgR);
      instance.SetDeleteArray(&deleteArray_vectorlEphyscLcLDiBosonlEphyscLcLParticlecOphyscLcLParticlegRsPgR);
      instance.SetDestructor(&destruct_vectorlEphyscLcLDiBosonlEphyscLcLParticlecOphyscLcLParticlegRsPgR);
      instance.AdoptCollectionProxyInfo(TCollectionProxyInfo::Generate(TCollectionProxyInfo::Pushback< vector<phys::DiBoson<phys::Particle,phys::Particle> > >()));
      return &instance;
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal((const vector<phys::DiBoson<phys::Particle,phys::Particle> >*)0x0); R__UseDummy(_R__UNIQUE_DICT_(Init));

   // Dictionary for non-ClassDef classes
   static TClass *vectorlEphyscLcLDiBosonlEphyscLcLParticlecOphyscLcLParticlegRsPgR_Dictionary() {
      TClass* theClass =::ROOT::GenerateInitInstanceLocal((const vector<phys::DiBoson<phys::Particle,phys::Particle> >*)0x0)->GetClass();
      vectorlEphyscLcLDiBosonlEphyscLcLParticlecOphyscLcLParticlegRsPgR_TClassManip(theClass);
   return theClass;
   }

   static void vectorlEphyscLcLDiBosonlEphyscLcLParticlecOphyscLcLParticlegRsPgR_TClassManip(TClass* ){
   }

} // end of namespace ROOT

namespace ROOT {
   // Wrappers around operator new
   static void *new_vectorlEphyscLcLDiBosonlEphyscLcLParticlecOphyscLcLParticlegRsPgR(void *p) {
      return  p ? ::new((::ROOT::Internal::TOperatorNewHelper*)p) vector<phys::DiBoson<phys::Particle,phys::Particle> > : new vector<phys::DiBoson<phys::Particle,phys::Particle> >;
   }
   static void *newArray_vectorlEphyscLcLDiBosonlEphyscLcLParticlecOphyscLcLParticlegRsPgR(Long_t nElements, void *p) {
      return p ? ::new((::ROOT::Internal::TOperatorNewHelper*)p) vector<phys::DiBoson<phys::Particle,phys::Particle> >[nElements] : new vector<phys::DiBoson<phys::Particle,phys::Particle> >[nElements];
   }
   // Wrapper around operator delete
   static void delete_vectorlEphyscLcLDiBosonlEphyscLcLParticlecOphyscLcLParticlegRsPgR(void *p) {
      delete ((vector<phys::DiBoson<phys::Particle,phys::Particle> >*)p);
   }
   static void deleteArray_vectorlEphyscLcLDiBosonlEphyscLcLParticlecOphyscLcLParticlegRsPgR(void *p) {
      delete [] ((vector<phys::DiBoson<phys::Particle,phys::Particle> >*)p);
   }
   static void destruct_vectorlEphyscLcLDiBosonlEphyscLcLParticlecOphyscLcLParticlegRsPgR(void *p) {
      typedef vector<phys::DiBoson<phys::Particle,phys::Particle> > current_t;
      ((current_t*)p)->~current_t();
   }
} // end of namespace ROOT for class vector<phys::DiBoson<phys::Particle,phys::Particle> >

namespace ROOT {
   static TClass *vectorlEphyscLcLDiBosonlEphyscLcLLeptoncOphyscLcLLeptongRsPgR_Dictionary();
   static void vectorlEphyscLcLDiBosonlEphyscLcLLeptoncOphyscLcLLeptongRsPgR_TClassManip(TClass*);
   static void *new_vectorlEphyscLcLDiBosonlEphyscLcLLeptoncOphyscLcLLeptongRsPgR(void *p = 0);
   static void *newArray_vectorlEphyscLcLDiBosonlEphyscLcLLeptoncOphyscLcLLeptongRsPgR(Long_t size, void *p);
   static void delete_vectorlEphyscLcLDiBosonlEphyscLcLLeptoncOphyscLcLLeptongRsPgR(void *p);
   static void deleteArray_vectorlEphyscLcLDiBosonlEphyscLcLLeptoncOphyscLcLLeptongRsPgR(void *p);
   static void destruct_vectorlEphyscLcLDiBosonlEphyscLcLLeptoncOphyscLcLLeptongRsPgR(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const vector<phys::DiBoson<phys::Lepton,phys::Lepton> >*)
   {
      vector<phys::DiBoson<phys::Lepton,phys::Lepton> > *ptr = 0;
      static ::TVirtualIsAProxy* isa_proxy = new ::TIsAProxy(typeid(vector<phys::DiBoson<phys::Lepton,phys::Lepton> >));
      static ::ROOT::TGenericClassInfo 
         instance("vector<phys::DiBoson<phys::Lepton,phys::Lepton> >", -2, "vector", 214,
                  typeid(vector<phys::DiBoson<phys::Lepton,phys::Lepton> >), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &vectorlEphyscLcLDiBosonlEphyscLcLLeptoncOphyscLcLLeptongRsPgR_Dictionary, isa_proxy, 4,
                  sizeof(vector<phys::DiBoson<phys::Lepton,phys::Lepton> >) );
      instance.SetNew(&new_vectorlEphyscLcLDiBosonlEphyscLcLLeptoncOphyscLcLLeptongRsPgR);
      instance.SetNewArray(&newArray_vectorlEphyscLcLDiBosonlEphyscLcLLeptoncOphyscLcLLeptongRsPgR);
      instance.SetDelete(&delete_vectorlEphyscLcLDiBosonlEphyscLcLLeptoncOphyscLcLLeptongRsPgR);
      instance.SetDeleteArray(&deleteArray_vectorlEphyscLcLDiBosonlEphyscLcLLeptoncOphyscLcLLeptongRsPgR);
      instance.SetDestructor(&destruct_vectorlEphyscLcLDiBosonlEphyscLcLLeptoncOphyscLcLLeptongRsPgR);
      instance.AdoptCollectionProxyInfo(TCollectionProxyInfo::Generate(TCollectionProxyInfo::Pushback< vector<phys::DiBoson<phys::Lepton,phys::Lepton> > >()));
      return &instance;
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal((const vector<phys::DiBoson<phys::Lepton,phys::Lepton> >*)0x0); R__UseDummy(_R__UNIQUE_DICT_(Init));

   // Dictionary for non-ClassDef classes
   static TClass *vectorlEphyscLcLDiBosonlEphyscLcLLeptoncOphyscLcLLeptongRsPgR_Dictionary() {
      TClass* theClass =::ROOT::GenerateInitInstanceLocal((const vector<phys::DiBoson<phys::Lepton,phys::Lepton> >*)0x0)->GetClass();
      vectorlEphyscLcLDiBosonlEphyscLcLLeptoncOphyscLcLLeptongRsPgR_TClassManip(theClass);
   return theClass;
   }

   static void vectorlEphyscLcLDiBosonlEphyscLcLLeptoncOphyscLcLLeptongRsPgR_TClassManip(TClass* ){
   }

} // end of namespace ROOT

namespace ROOT {
   // Wrappers around operator new
   static void *new_vectorlEphyscLcLDiBosonlEphyscLcLLeptoncOphyscLcLLeptongRsPgR(void *p) {
      return  p ? ::new((::ROOT::Internal::TOperatorNewHelper*)p) vector<phys::DiBoson<phys::Lepton,phys::Lepton> > : new vector<phys::DiBoson<phys::Lepton,phys::Lepton> >;
   }
   static void *newArray_vectorlEphyscLcLDiBosonlEphyscLcLLeptoncOphyscLcLLeptongRsPgR(Long_t nElements, void *p) {
      return p ? ::new((::ROOT::Internal::TOperatorNewHelper*)p) vector<phys::DiBoson<phys::Lepton,phys::Lepton> >[nElements] : new vector<phys::DiBoson<phys::Lepton,phys::Lepton> >[nElements];
   }
   // Wrapper around operator delete
   static void delete_vectorlEphyscLcLDiBosonlEphyscLcLLeptoncOphyscLcLLeptongRsPgR(void *p) {
      delete ((vector<phys::DiBoson<phys::Lepton,phys::Lepton> >*)p);
   }
   static void deleteArray_vectorlEphyscLcLDiBosonlEphyscLcLLeptoncOphyscLcLLeptongRsPgR(void *p) {
      delete [] ((vector<phys::DiBoson<phys::Lepton,phys::Lepton> >*)p);
   }
   static void destruct_vectorlEphyscLcLDiBosonlEphyscLcLLeptoncOphyscLcLLeptongRsPgR(void *p) {
      typedef vector<phys::DiBoson<phys::Lepton,phys::Lepton> > current_t;
      ((current_t*)p)->~current_t();
   }
} // end of namespace ROOT for class vector<phys::DiBoson<phys::Lepton,phys::Lepton> >

namespace ROOT {
   static TClass *vectorlEphyscLcLDiBosonlEphyscLcLLeptoncOphyscLcLElectrongRsPgR_Dictionary();
   static void vectorlEphyscLcLDiBosonlEphyscLcLLeptoncOphyscLcLElectrongRsPgR_TClassManip(TClass*);
   static void *new_vectorlEphyscLcLDiBosonlEphyscLcLLeptoncOphyscLcLElectrongRsPgR(void *p = 0);
   static void *newArray_vectorlEphyscLcLDiBosonlEphyscLcLLeptoncOphyscLcLElectrongRsPgR(Long_t size, void *p);
   static void delete_vectorlEphyscLcLDiBosonlEphyscLcLLeptoncOphyscLcLElectrongRsPgR(void *p);
   static void deleteArray_vectorlEphyscLcLDiBosonlEphyscLcLLeptoncOphyscLcLElectrongRsPgR(void *p);
   static void destruct_vectorlEphyscLcLDiBosonlEphyscLcLLeptoncOphyscLcLElectrongRsPgR(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const vector<phys::DiBoson<phys::Lepton,phys::Electron> >*)
   {
      vector<phys::DiBoson<phys::Lepton,phys::Electron> > *ptr = 0;
      static ::TVirtualIsAProxy* isa_proxy = new ::TIsAProxy(typeid(vector<phys::DiBoson<phys::Lepton,phys::Electron> >));
      static ::ROOT::TGenericClassInfo 
         instance("vector<phys::DiBoson<phys::Lepton,phys::Electron> >", -2, "vector", 214,
                  typeid(vector<phys::DiBoson<phys::Lepton,phys::Electron> >), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &vectorlEphyscLcLDiBosonlEphyscLcLLeptoncOphyscLcLElectrongRsPgR_Dictionary, isa_proxy, 4,
                  sizeof(vector<phys::DiBoson<phys::Lepton,phys::Electron> >) );
      instance.SetNew(&new_vectorlEphyscLcLDiBosonlEphyscLcLLeptoncOphyscLcLElectrongRsPgR);
      instance.SetNewArray(&newArray_vectorlEphyscLcLDiBosonlEphyscLcLLeptoncOphyscLcLElectrongRsPgR);
      instance.SetDelete(&delete_vectorlEphyscLcLDiBosonlEphyscLcLLeptoncOphyscLcLElectrongRsPgR);
      instance.SetDeleteArray(&deleteArray_vectorlEphyscLcLDiBosonlEphyscLcLLeptoncOphyscLcLElectrongRsPgR);
      instance.SetDestructor(&destruct_vectorlEphyscLcLDiBosonlEphyscLcLLeptoncOphyscLcLElectrongRsPgR);
      instance.AdoptCollectionProxyInfo(TCollectionProxyInfo::Generate(TCollectionProxyInfo::Pushback< vector<phys::DiBoson<phys::Lepton,phys::Electron> > >()));
      return &instance;
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal((const vector<phys::DiBoson<phys::Lepton,phys::Electron> >*)0x0); R__UseDummy(_R__UNIQUE_DICT_(Init));

   // Dictionary for non-ClassDef classes
   static TClass *vectorlEphyscLcLDiBosonlEphyscLcLLeptoncOphyscLcLElectrongRsPgR_Dictionary() {
      TClass* theClass =::ROOT::GenerateInitInstanceLocal((const vector<phys::DiBoson<phys::Lepton,phys::Electron> >*)0x0)->GetClass();
      vectorlEphyscLcLDiBosonlEphyscLcLLeptoncOphyscLcLElectrongRsPgR_TClassManip(theClass);
   return theClass;
   }

   static void vectorlEphyscLcLDiBosonlEphyscLcLLeptoncOphyscLcLElectrongRsPgR_TClassManip(TClass* ){
   }

} // end of namespace ROOT

namespace ROOT {
   // Wrappers around operator new
   static void *new_vectorlEphyscLcLDiBosonlEphyscLcLLeptoncOphyscLcLElectrongRsPgR(void *p) {
      return  p ? ::new((::ROOT::Internal::TOperatorNewHelper*)p) vector<phys::DiBoson<phys::Lepton,phys::Electron> > : new vector<phys::DiBoson<phys::Lepton,phys::Electron> >;
   }
   static void *newArray_vectorlEphyscLcLDiBosonlEphyscLcLLeptoncOphyscLcLElectrongRsPgR(Long_t nElements, void *p) {
      return p ? ::new((::ROOT::Internal::TOperatorNewHelper*)p) vector<phys::DiBoson<phys::Lepton,phys::Electron> >[nElements] : new vector<phys::DiBoson<phys::Lepton,phys::Electron> >[nElements];
   }
   // Wrapper around operator delete
   static void delete_vectorlEphyscLcLDiBosonlEphyscLcLLeptoncOphyscLcLElectrongRsPgR(void *p) {
      delete ((vector<phys::DiBoson<phys::Lepton,phys::Electron> >*)p);
   }
   static void deleteArray_vectorlEphyscLcLDiBosonlEphyscLcLLeptoncOphyscLcLElectrongRsPgR(void *p) {
      delete [] ((vector<phys::DiBoson<phys::Lepton,phys::Electron> >*)p);
   }
   static void destruct_vectorlEphyscLcLDiBosonlEphyscLcLLeptoncOphyscLcLElectrongRsPgR(void *p) {
      typedef vector<phys::DiBoson<phys::Lepton,phys::Electron> > current_t;
      ((current_t*)p)->~current_t();
   }
} // end of namespace ROOT for class vector<phys::DiBoson<phys::Lepton,phys::Electron> >

namespace ROOT {
   static TClass *vectorlEphyscLcLDiBosonlEphyscLcLElectroncOphyscLcLLeptongRsPgR_Dictionary();
   static void vectorlEphyscLcLDiBosonlEphyscLcLElectroncOphyscLcLLeptongRsPgR_TClassManip(TClass*);
   static void *new_vectorlEphyscLcLDiBosonlEphyscLcLElectroncOphyscLcLLeptongRsPgR(void *p = 0);
   static void *newArray_vectorlEphyscLcLDiBosonlEphyscLcLElectroncOphyscLcLLeptongRsPgR(Long_t size, void *p);
   static void delete_vectorlEphyscLcLDiBosonlEphyscLcLElectroncOphyscLcLLeptongRsPgR(void *p);
   static void deleteArray_vectorlEphyscLcLDiBosonlEphyscLcLElectroncOphyscLcLLeptongRsPgR(void *p);
   static void destruct_vectorlEphyscLcLDiBosonlEphyscLcLElectroncOphyscLcLLeptongRsPgR(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const vector<phys::DiBoson<phys::Electron,phys::Lepton> >*)
   {
      vector<phys::DiBoson<phys::Electron,phys::Lepton> > *ptr = 0;
      static ::TVirtualIsAProxy* isa_proxy = new ::TIsAProxy(typeid(vector<phys::DiBoson<phys::Electron,phys::Lepton> >));
      static ::ROOT::TGenericClassInfo 
         instance("vector<phys::DiBoson<phys::Electron,phys::Lepton> >", -2, "vector", 214,
                  typeid(vector<phys::DiBoson<phys::Electron,phys::Lepton> >), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &vectorlEphyscLcLDiBosonlEphyscLcLElectroncOphyscLcLLeptongRsPgR_Dictionary, isa_proxy, 4,
                  sizeof(vector<phys::DiBoson<phys::Electron,phys::Lepton> >) );
      instance.SetNew(&new_vectorlEphyscLcLDiBosonlEphyscLcLElectroncOphyscLcLLeptongRsPgR);
      instance.SetNewArray(&newArray_vectorlEphyscLcLDiBosonlEphyscLcLElectroncOphyscLcLLeptongRsPgR);
      instance.SetDelete(&delete_vectorlEphyscLcLDiBosonlEphyscLcLElectroncOphyscLcLLeptongRsPgR);
      instance.SetDeleteArray(&deleteArray_vectorlEphyscLcLDiBosonlEphyscLcLElectroncOphyscLcLLeptongRsPgR);
      instance.SetDestructor(&destruct_vectorlEphyscLcLDiBosonlEphyscLcLElectroncOphyscLcLLeptongRsPgR);
      instance.AdoptCollectionProxyInfo(TCollectionProxyInfo::Generate(TCollectionProxyInfo::Pushback< vector<phys::DiBoson<phys::Electron,phys::Lepton> > >()));
      return &instance;
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal((const vector<phys::DiBoson<phys::Electron,phys::Lepton> >*)0x0); R__UseDummy(_R__UNIQUE_DICT_(Init));

   // Dictionary for non-ClassDef classes
   static TClass *vectorlEphyscLcLDiBosonlEphyscLcLElectroncOphyscLcLLeptongRsPgR_Dictionary() {
      TClass* theClass =::ROOT::GenerateInitInstanceLocal((const vector<phys::DiBoson<phys::Electron,phys::Lepton> >*)0x0)->GetClass();
      vectorlEphyscLcLDiBosonlEphyscLcLElectroncOphyscLcLLeptongRsPgR_TClassManip(theClass);
   return theClass;
   }

   static void vectorlEphyscLcLDiBosonlEphyscLcLElectroncOphyscLcLLeptongRsPgR_TClassManip(TClass* ){
   }

} // end of namespace ROOT

namespace ROOT {
   // Wrappers around operator new
   static void *new_vectorlEphyscLcLDiBosonlEphyscLcLElectroncOphyscLcLLeptongRsPgR(void *p) {
      return  p ? ::new((::ROOT::Internal::TOperatorNewHelper*)p) vector<phys::DiBoson<phys::Electron,phys::Lepton> > : new vector<phys::DiBoson<phys::Electron,phys::Lepton> >;
   }
   static void *newArray_vectorlEphyscLcLDiBosonlEphyscLcLElectroncOphyscLcLLeptongRsPgR(Long_t nElements, void *p) {
      return p ? ::new((::ROOT::Internal::TOperatorNewHelper*)p) vector<phys::DiBoson<phys::Electron,phys::Lepton> >[nElements] : new vector<phys::DiBoson<phys::Electron,phys::Lepton> >[nElements];
   }
   // Wrapper around operator delete
   static void delete_vectorlEphyscLcLDiBosonlEphyscLcLElectroncOphyscLcLLeptongRsPgR(void *p) {
      delete ((vector<phys::DiBoson<phys::Electron,phys::Lepton> >*)p);
   }
   static void deleteArray_vectorlEphyscLcLDiBosonlEphyscLcLElectroncOphyscLcLLeptongRsPgR(void *p) {
      delete [] ((vector<phys::DiBoson<phys::Electron,phys::Lepton> >*)p);
   }
   static void destruct_vectorlEphyscLcLDiBosonlEphyscLcLElectroncOphyscLcLLeptongRsPgR(void *p) {
      typedef vector<phys::DiBoson<phys::Electron,phys::Lepton> > current_t;
      ((current_t*)p)->~current_t();
   }
} // end of namespace ROOT for class vector<phys::DiBoson<phys::Electron,phys::Lepton> >

namespace ROOT {
   static TClass *vectorlEphyscLcLDiBosonlEphyscLcLElectroncOphyscLcLElectrongRsPgR_Dictionary();
   static void vectorlEphyscLcLDiBosonlEphyscLcLElectroncOphyscLcLElectrongRsPgR_TClassManip(TClass*);
   static void *new_vectorlEphyscLcLDiBosonlEphyscLcLElectroncOphyscLcLElectrongRsPgR(void *p = 0);
   static void *newArray_vectorlEphyscLcLDiBosonlEphyscLcLElectroncOphyscLcLElectrongRsPgR(Long_t size, void *p);
   static void delete_vectorlEphyscLcLDiBosonlEphyscLcLElectroncOphyscLcLElectrongRsPgR(void *p);
   static void deleteArray_vectorlEphyscLcLDiBosonlEphyscLcLElectroncOphyscLcLElectrongRsPgR(void *p);
   static void destruct_vectorlEphyscLcLDiBosonlEphyscLcLElectroncOphyscLcLElectrongRsPgR(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const vector<phys::DiBoson<phys::Electron,phys::Electron> >*)
   {
      vector<phys::DiBoson<phys::Electron,phys::Electron> > *ptr = 0;
      static ::TVirtualIsAProxy* isa_proxy = new ::TIsAProxy(typeid(vector<phys::DiBoson<phys::Electron,phys::Electron> >));
      static ::ROOT::TGenericClassInfo 
         instance("vector<phys::DiBoson<phys::Electron,phys::Electron> >", -2, "vector", 214,
                  typeid(vector<phys::DiBoson<phys::Electron,phys::Electron> >), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &vectorlEphyscLcLDiBosonlEphyscLcLElectroncOphyscLcLElectrongRsPgR_Dictionary, isa_proxy, 4,
                  sizeof(vector<phys::DiBoson<phys::Electron,phys::Electron> >) );
      instance.SetNew(&new_vectorlEphyscLcLDiBosonlEphyscLcLElectroncOphyscLcLElectrongRsPgR);
      instance.SetNewArray(&newArray_vectorlEphyscLcLDiBosonlEphyscLcLElectroncOphyscLcLElectrongRsPgR);
      instance.SetDelete(&delete_vectorlEphyscLcLDiBosonlEphyscLcLElectroncOphyscLcLElectrongRsPgR);
      instance.SetDeleteArray(&deleteArray_vectorlEphyscLcLDiBosonlEphyscLcLElectroncOphyscLcLElectrongRsPgR);
      instance.SetDestructor(&destruct_vectorlEphyscLcLDiBosonlEphyscLcLElectroncOphyscLcLElectrongRsPgR);
      instance.AdoptCollectionProxyInfo(TCollectionProxyInfo::Generate(TCollectionProxyInfo::Pushback< vector<phys::DiBoson<phys::Electron,phys::Electron> > >()));
      return &instance;
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal((const vector<phys::DiBoson<phys::Electron,phys::Electron> >*)0x0); R__UseDummy(_R__UNIQUE_DICT_(Init));

   // Dictionary for non-ClassDef classes
   static TClass *vectorlEphyscLcLDiBosonlEphyscLcLElectroncOphyscLcLElectrongRsPgR_Dictionary() {
      TClass* theClass =::ROOT::GenerateInitInstanceLocal((const vector<phys::DiBoson<phys::Electron,phys::Electron> >*)0x0)->GetClass();
      vectorlEphyscLcLDiBosonlEphyscLcLElectroncOphyscLcLElectrongRsPgR_TClassManip(theClass);
   return theClass;
   }

   static void vectorlEphyscLcLDiBosonlEphyscLcLElectroncOphyscLcLElectrongRsPgR_TClassManip(TClass* ){
   }

} // end of namespace ROOT

namespace ROOT {
   // Wrappers around operator new
   static void *new_vectorlEphyscLcLDiBosonlEphyscLcLElectroncOphyscLcLElectrongRsPgR(void *p) {
      return  p ? ::new((::ROOT::Internal::TOperatorNewHelper*)p) vector<phys::DiBoson<phys::Electron,phys::Electron> > : new vector<phys::DiBoson<phys::Electron,phys::Electron> >;
   }
   static void *newArray_vectorlEphyscLcLDiBosonlEphyscLcLElectroncOphyscLcLElectrongRsPgR(Long_t nElements, void *p) {
      return p ? ::new((::ROOT::Internal::TOperatorNewHelper*)p) vector<phys::DiBoson<phys::Electron,phys::Electron> >[nElements] : new vector<phys::DiBoson<phys::Electron,phys::Electron> >[nElements];
   }
   // Wrapper around operator delete
   static void delete_vectorlEphyscLcLDiBosonlEphyscLcLElectroncOphyscLcLElectrongRsPgR(void *p) {
      delete ((vector<phys::DiBoson<phys::Electron,phys::Electron> >*)p);
   }
   static void deleteArray_vectorlEphyscLcLDiBosonlEphyscLcLElectroncOphyscLcLElectrongRsPgR(void *p) {
      delete [] ((vector<phys::DiBoson<phys::Electron,phys::Electron> >*)p);
   }
   static void destruct_vectorlEphyscLcLDiBosonlEphyscLcLElectroncOphyscLcLElectrongRsPgR(void *p) {
      typedef vector<phys::DiBoson<phys::Electron,phys::Electron> > current_t;
      ((current_t*)p)->~current_t();
   }
} // end of namespace ROOT for class vector<phys::DiBoson<phys::Electron,phys::Electron> >

namespace ROOT {
   static TClass *vectorlEphyscLcLBosonlEphyscLcLParticlegRsPgR_Dictionary();
   static void vectorlEphyscLcLBosonlEphyscLcLParticlegRsPgR_TClassManip(TClass*);
   static void *new_vectorlEphyscLcLBosonlEphyscLcLParticlegRsPgR(void *p = 0);
   static void *newArray_vectorlEphyscLcLBosonlEphyscLcLParticlegRsPgR(Long_t size, void *p);
   static void delete_vectorlEphyscLcLBosonlEphyscLcLParticlegRsPgR(void *p);
   static void deleteArray_vectorlEphyscLcLBosonlEphyscLcLParticlegRsPgR(void *p);
   static void destruct_vectorlEphyscLcLBosonlEphyscLcLParticlegRsPgR(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const vector<phys::Boson<phys::Particle> >*)
   {
      vector<phys::Boson<phys::Particle> > *ptr = 0;
      static ::TVirtualIsAProxy* isa_proxy = new ::TIsAProxy(typeid(vector<phys::Boson<phys::Particle> >));
      static ::ROOT::TGenericClassInfo 
         instance("vector<phys::Boson<phys::Particle> >", -2, "vector", 214,
                  typeid(vector<phys::Boson<phys::Particle> >), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &vectorlEphyscLcLBosonlEphyscLcLParticlegRsPgR_Dictionary, isa_proxy, 4,
                  sizeof(vector<phys::Boson<phys::Particle> >) );
      instance.SetNew(&new_vectorlEphyscLcLBosonlEphyscLcLParticlegRsPgR);
      instance.SetNewArray(&newArray_vectorlEphyscLcLBosonlEphyscLcLParticlegRsPgR);
      instance.SetDelete(&delete_vectorlEphyscLcLBosonlEphyscLcLParticlegRsPgR);
      instance.SetDeleteArray(&deleteArray_vectorlEphyscLcLBosonlEphyscLcLParticlegRsPgR);
      instance.SetDestructor(&destruct_vectorlEphyscLcLBosonlEphyscLcLParticlegRsPgR);
      instance.AdoptCollectionProxyInfo(TCollectionProxyInfo::Generate(TCollectionProxyInfo::Pushback< vector<phys::Boson<phys::Particle> > >()));
      return &instance;
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal((const vector<phys::Boson<phys::Particle> >*)0x0); R__UseDummy(_R__UNIQUE_DICT_(Init));

   // Dictionary for non-ClassDef classes
   static TClass *vectorlEphyscLcLBosonlEphyscLcLParticlegRsPgR_Dictionary() {
      TClass* theClass =::ROOT::GenerateInitInstanceLocal((const vector<phys::Boson<phys::Particle> >*)0x0)->GetClass();
      vectorlEphyscLcLBosonlEphyscLcLParticlegRsPgR_TClassManip(theClass);
   return theClass;
   }

   static void vectorlEphyscLcLBosonlEphyscLcLParticlegRsPgR_TClassManip(TClass* ){
   }

} // end of namespace ROOT

namespace ROOT {
   // Wrappers around operator new
   static void *new_vectorlEphyscLcLBosonlEphyscLcLParticlegRsPgR(void *p) {
      return  p ? ::new((::ROOT::Internal::TOperatorNewHelper*)p) vector<phys::Boson<phys::Particle> > : new vector<phys::Boson<phys::Particle> >;
   }
   static void *newArray_vectorlEphyscLcLBosonlEphyscLcLParticlegRsPgR(Long_t nElements, void *p) {
      return p ? ::new((::ROOT::Internal::TOperatorNewHelper*)p) vector<phys::Boson<phys::Particle> >[nElements] : new vector<phys::Boson<phys::Particle> >[nElements];
   }
   // Wrapper around operator delete
   static void delete_vectorlEphyscLcLBosonlEphyscLcLParticlegRsPgR(void *p) {
      delete ((vector<phys::Boson<phys::Particle> >*)p);
   }
   static void deleteArray_vectorlEphyscLcLBosonlEphyscLcLParticlegRsPgR(void *p) {
      delete [] ((vector<phys::Boson<phys::Particle> >*)p);
   }
   static void destruct_vectorlEphyscLcLBosonlEphyscLcLParticlegRsPgR(void *p) {
      typedef vector<phys::Boson<phys::Particle> > current_t;
      ((current_t*)p)->~current_t();
   }
} // end of namespace ROOT for class vector<phys::Boson<phys::Particle> >

namespace ROOT {
   static TClass *vectorlEphyscLcLBosonlEphyscLcLLeptongRsPgR_Dictionary();
   static void vectorlEphyscLcLBosonlEphyscLcLLeptongRsPgR_TClassManip(TClass*);
   static void *new_vectorlEphyscLcLBosonlEphyscLcLLeptongRsPgR(void *p = 0);
   static void *newArray_vectorlEphyscLcLBosonlEphyscLcLLeptongRsPgR(Long_t size, void *p);
   static void delete_vectorlEphyscLcLBosonlEphyscLcLLeptongRsPgR(void *p);
   static void deleteArray_vectorlEphyscLcLBosonlEphyscLcLLeptongRsPgR(void *p);
   static void destruct_vectorlEphyscLcLBosonlEphyscLcLLeptongRsPgR(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const vector<phys::Boson<phys::Lepton> >*)
   {
      vector<phys::Boson<phys::Lepton> > *ptr = 0;
      static ::TVirtualIsAProxy* isa_proxy = new ::TIsAProxy(typeid(vector<phys::Boson<phys::Lepton> >));
      static ::ROOT::TGenericClassInfo 
         instance("vector<phys::Boson<phys::Lepton> >", -2, "vector", 214,
                  typeid(vector<phys::Boson<phys::Lepton> >), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &vectorlEphyscLcLBosonlEphyscLcLLeptongRsPgR_Dictionary, isa_proxy, 4,
                  sizeof(vector<phys::Boson<phys::Lepton> >) );
      instance.SetNew(&new_vectorlEphyscLcLBosonlEphyscLcLLeptongRsPgR);
      instance.SetNewArray(&newArray_vectorlEphyscLcLBosonlEphyscLcLLeptongRsPgR);
      instance.SetDelete(&delete_vectorlEphyscLcLBosonlEphyscLcLLeptongRsPgR);
      instance.SetDeleteArray(&deleteArray_vectorlEphyscLcLBosonlEphyscLcLLeptongRsPgR);
      instance.SetDestructor(&destruct_vectorlEphyscLcLBosonlEphyscLcLLeptongRsPgR);
      instance.AdoptCollectionProxyInfo(TCollectionProxyInfo::Generate(TCollectionProxyInfo::Pushback< vector<phys::Boson<phys::Lepton> > >()));
      return &instance;
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal((const vector<phys::Boson<phys::Lepton> >*)0x0); R__UseDummy(_R__UNIQUE_DICT_(Init));

   // Dictionary for non-ClassDef classes
   static TClass *vectorlEphyscLcLBosonlEphyscLcLLeptongRsPgR_Dictionary() {
      TClass* theClass =::ROOT::GenerateInitInstanceLocal((const vector<phys::Boson<phys::Lepton> >*)0x0)->GetClass();
      vectorlEphyscLcLBosonlEphyscLcLLeptongRsPgR_TClassManip(theClass);
   return theClass;
   }

   static void vectorlEphyscLcLBosonlEphyscLcLLeptongRsPgR_TClassManip(TClass* ){
   }

} // end of namespace ROOT

namespace ROOT {
   // Wrappers around operator new
   static void *new_vectorlEphyscLcLBosonlEphyscLcLLeptongRsPgR(void *p) {
      return  p ? ::new((::ROOT::Internal::TOperatorNewHelper*)p) vector<phys::Boson<phys::Lepton> > : new vector<phys::Boson<phys::Lepton> >;
   }
   static void *newArray_vectorlEphyscLcLBosonlEphyscLcLLeptongRsPgR(Long_t nElements, void *p) {
      return p ? ::new((::ROOT::Internal::TOperatorNewHelper*)p) vector<phys::Boson<phys::Lepton> >[nElements] : new vector<phys::Boson<phys::Lepton> >[nElements];
   }
   // Wrapper around operator delete
   static void delete_vectorlEphyscLcLBosonlEphyscLcLLeptongRsPgR(void *p) {
      delete ((vector<phys::Boson<phys::Lepton> >*)p);
   }
   static void deleteArray_vectorlEphyscLcLBosonlEphyscLcLLeptongRsPgR(void *p) {
      delete [] ((vector<phys::Boson<phys::Lepton> >*)p);
   }
   static void destruct_vectorlEphyscLcLBosonlEphyscLcLLeptongRsPgR(void *p) {
      typedef vector<phys::Boson<phys::Lepton> > current_t;
      ((current_t*)p)->~current_t();
   }
} // end of namespace ROOT for class vector<phys::Boson<phys::Lepton> >

namespace ROOT {
   static TClass *vectorlEphyscLcLBosonlEphyscLcLJetgRsPgR_Dictionary();
   static void vectorlEphyscLcLBosonlEphyscLcLJetgRsPgR_TClassManip(TClass*);
   static void *new_vectorlEphyscLcLBosonlEphyscLcLJetgRsPgR(void *p = 0);
   static void *newArray_vectorlEphyscLcLBosonlEphyscLcLJetgRsPgR(Long_t size, void *p);
   static void delete_vectorlEphyscLcLBosonlEphyscLcLJetgRsPgR(void *p);
   static void deleteArray_vectorlEphyscLcLBosonlEphyscLcLJetgRsPgR(void *p);
   static void destruct_vectorlEphyscLcLBosonlEphyscLcLJetgRsPgR(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const vector<phys::Boson<phys::Jet> >*)
   {
      vector<phys::Boson<phys::Jet> > *ptr = 0;
      static ::TVirtualIsAProxy* isa_proxy = new ::TIsAProxy(typeid(vector<phys::Boson<phys::Jet> >));
      static ::ROOT::TGenericClassInfo 
         instance("vector<phys::Boson<phys::Jet> >", -2, "vector", 214,
                  typeid(vector<phys::Boson<phys::Jet> >), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &vectorlEphyscLcLBosonlEphyscLcLJetgRsPgR_Dictionary, isa_proxy, 4,
                  sizeof(vector<phys::Boson<phys::Jet> >) );
      instance.SetNew(&new_vectorlEphyscLcLBosonlEphyscLcLJetgRsPgR);
      instance.SetNewArray(&newArray_vectorlEphyscLcLBosonlEphyscLcLJetgRsPgR);
      instance.SetDelete(&delete_vectorlEphyscLcLBosonlEphyscLcLJetgRsPgR);
      instance.SetDeleteArray(&deleteArray_vectorlEphyscLcLBosonlEphyscLcLJetgRsPgR);
      instance.SetDestructor(&destruct_vectorlEphyscLcLBosonlEphyscLcLJetgRsPgR);
      instance.AdoptCollectionProxyInfo(TCollectionProxyInfo::Generate(TCollectionProxyInfo::Pushback< vector<phys::Boson<phys::Jet> > >()));
      return &instance;
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal((const vector<phys::Boson<phys::Jet> >*)0x0); R__UseDummy(_R__UNIQUE_DICT_(Init));

   // Dictionary for non-ClassDef classes
   static TClass *vectorlEphyscLcLBosonlEphyscLcLJetgRsPgR_Dictionary() {
      TClass* theClass =::ROOT::GenerateInitInstanceLocal((const vector<phys::Boson<phys::Jet> >*)0x0)->GetClass();
      vectorlEphyscLcLBosonlEphyscLcLJetgRsPgR_TClassManip(theClass);
   return theClass;
   }

   static void vectorlEphyscLcLBosonlEphyscLcLJetgRsPgR_TClassManip(TClass* ){
   }

} // end of namespace ROOT

namespace ROOT {
   // Wrappers around operator new
   static void *new_vectorlEphyscLcLBosonlEphyscLcLJetgRsPgR(void *p) {
      return  p ? ::new((::ROOT::Internal::TOperatorNewHelper*)p) vector<phys::Boson<phys::Jet> > : new vector<phys::Boson<phys::Jet> >;
   }
   static void *newArray_vectorlEphyscLcLBosonlEphyscLcLJetgRsPgR(Long_t nElements, void *p) {
      return p ? ::new((::ROOT::Internal::TOperatorNewHelper*)p) vector<phys::Boson<phys::Jet> >[nElements] : new vector<phys::Boson<phys::Jet> >[nElements];
   }
   // Wrapper around operator delete
   static void delete_vectorlEphyscLcLBosonlEphyscLcLJetgRsPgR(void *p) {
      delete ((vector<phys::Boson<phys::Jet> >*)p);
   }
   static void deleteArray_vectorlEphyscLcLBosonlEphyscLcLJetgRsPgR(void *p) {
      delete [] ((vector<phys::Boson<phys::Jet> >*)p);
   }
   static void destruct_vectorlEphyscLcLBosonlEphyscLcLJetgRsPgR(void *p) {
      typedef vector<phys::Boson<phys::Jet> > current_t;
      ((current_t*)p)->~current_t();
   }
} // end of namespace ROOT for class vector<phys::Boson<phys::Jet> >

namespace ROOT {
   static TClass *vectorlEphyscLcLBosonlEphyscLcLElectrongRsPgR_Dictionary();
   static void vectorlEphyscLcLBosonlEphyscLcLElectrongRsPgR_TClassManip(TClass*);
   static void *new_vectorlEphyscLcLBosonlEphyscLcLElectrongRsPgR(void *p = 0);
   static void *newArray_vectorlEphyscLcLBosonlEphyscLcLElectrongRsPgR(Long_t size, void *p);
   static void delete_vectorlEphyscLcLBosonlEphyscLcLElectrongRsPgR(void *p);
   static void deleteArray_vectorlEphyscLcLBosonlEphyscLcLElectrongRsPgR(void *p);
   static void destruct_vectorlEphyscLcLBosonlEphyscLcLElectrongRsPgR(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const vector<phys::Boson<phys::Electron> >*)
   {
      vector<phys::Boson<phys::Electron> > *ptr = 0;
      static ::TVirtualIsAProxy* isa_proxy = new ::TIsAProxy(typeid(vector<phys::Boson<phys::Electron> >));
      static ::ROOT::TGenericClassInfo 
         instance("vector<phys::Boson<phys::Electron> >", -2, "vector", 214,
                  typeid(vector<phys::Boson<phys::Electron> >), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &vectorlEphyscLcLBosonlEphyscLcLElectrongRsPgR_Dictionary, isa_proxy, 4,
                  sizeof(vector<phys::Boson<phys::Electron> >) );
      instance.SetNew(&new_vectorlEphyscLcLBosonlEphyscLcLElectrongRsPgR);
      instance.SetNewArray(&newArray_vectorlEphyscLcLBosonlEphyscLcLElectrongRsPgR);
      instance.SetDelete(&delete_vectorlEphyscLcLBosonlEphyscLcLElectrongRsPgR);
      instance.SetDeleteArray(&deleteArray_vectorlEphyscLcLBosonlEphyscLcLElectrongRsPgR);
      instance.SetDestructor(&destruct_vectorlEphyscLcLBosonlEphyscLcLElectrongRsPgR);
      instance.AdoptCollectionProxyInfo(TCollectionProxyInfo::Generate(TCollectionProxyInfo::Pushback< vector<phys::Boson<phys::Electron> > >()));
      return &instance;
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal((const vector<phys::Boson<phys::Electron> >*)0x0); R__UseDummy(_R__UNIQUE_DICT_(Init));

   // Dictionary for non-ClassDef classes
   static TClass *vectorlEphyscLcLBosonlEphyscLcLElectrongRsPgR_Dictionary() {
      TClass* theClass =::ROOT::GenerateInitInstanceLocal((const vector<phys::Boson<phys::Electron> >*)0x0)->GetClass();
      vectorlEphyscLcLBosonlEphyscLcLElectrongRsPgR_TClassManip(theClass);
   return theClass;
   }

   static void vectorlEphyscLcLBosonlEphyscLcLElectrongRsPgR_TClassManip(TClass* ){
   }

} // end of namespace ROOT

namespace ROOT {
   // Wrappers around operator new
   static void *new_vectorlEphyscLcLBosonlEphyscLcLElectrongRsPgR(void *p) {
      return  p ? ::new((::ROOT::Internal::TOperatorNewHelper*)p) vector<phys::Boson<phys::Electron> > : new vector<phys::Boson<phys::Electron> >;
   }
   static void *newArray_vectorlEphyscLcLBosonlEphyscLcLElectrongRsPgR(Long_t nElements, void *p) {
      return p ? ::new((::ROOT::Internal::TOperatorNewHelper*)p) vector<phys::Boson<phys::Electron> >[nElements] : new vector<phys::Boson<phys::Electron> >[nElements];
   }
   // Wrapper around operator delete
   static void delete_vectorlEphyscLcLBosonlEphyscLcLElectrongRsPgR(void *p) {
      delete ((vector<phys::Boson<phys::Electron> >*)p);
   }
   static void deleteArray_vectorlEphyscLcLBosonlEphyscLcLElectrongRsPgR(void *p) {
      delete [] ((vector<phys::Boson<phys::Electron> >*)p);
   }
   static void destruct_vectorlEphyscLcLBosonlEphyscLcLElectrongRsPgR(void *p) {
      typedef vector<phys::Boson<phys::Electron> > current_t;
      ((current_t*)p)->~current_t();
   }
} // end of namespace ROOT for class vector<phys::Boson<phys::Electron> >

namespace ROOT {
   static TClass *vectorlEpairlEphyscLcLBosonlEphyscLcLLeptongRcOphyscLcLLeptongRsPgR_Dictionary();
   static void vectorlEpairlEphyscLcLBosonlEphyscLcLLeptongRcOphyscLcLLeptongRsPgR_TClassManip(TClass*);
   static void *new_vectorlEpairlEphyscLcLBosonlEphyscLcLLeptongRcOphyscLcLLeptongRsPgR(void *p = 0);
   static void *newArray_vectorlEpairlEphyscLcLBosonlEphyscLcLLeptongRcOphyscLcLLeptongRsPgR(Long_t size, void *p);
   static void delete_vectorlEpairlEphyscLcLBosonlEphyscLcLLeptongRcOphyscLcLLeptongRsPgR(void *p);
   static void deleteArray_vectorlEpairlEphyscLcLBosonlEphyscLcLLeptongRcOphyscLcLLeptongRsPgR(void *p);
   static void destruct_vectorlEpairlEphyscLcLBosonlEphyscLcLLeptongRcOphyscLcLLeptongRsPgR(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const vector<pair<phys::Boson<phys::Lepton>,phys::Lepton> >*)
   {
      vector<pair<phys::Boson<phys::Lepton>,phys::Lepton> > *ptr = 0;
      static ::TVirtualIsAProxy* isa_proxy = new ::TIsAProxy(typeid(vector<pair<phys::Boson<phys::Lepton>,phys::Lepton> >));
      static ::ROOT::TGenericClassInfo 
         instance("vector<pair<phys::Boson<phys::Lepton>,phys::Lepton> >", -2, "vector", 214,
                  typeid(vector<pair<phys::Boson<phys::Lepton>,phys::Lepton> >), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &vectorlEpairlEphyscLcLBosonlEphyscLcLLeptongRcOphyscLcLLeptongRsPgR_Dictionary, isa_proxy, 4,
                  sizeof(vector<pair<phys::Boson<phys::Lepton>,phys::Lepton> >) );
      instance.SetNew(&new_vectorlEpairlEphyscLcLBosonlEphyscLcLLeptongRcOphyscLcLLeptongRsPgR);
      instance.SetNewArray(&newArray_vectorlEpairlEphyscLcLBosonlEphyscLcLLeptongRcOphyscLcLLeptongRsPgR);
      instance.SetDelete(&delete_vectorlEpairlEphyscLcLBosonlEphyscLcLLeptongRcOphyscLcLLeptongRsPgR);
      instance.SetDeleteArray(&deleteArray_vectorlEpairlEphyscLcLBosonlEphyscLcLLeptongRcOphyscLcLLeptongRsPgR);
      instance.SetDestructor(&destruct_vectorlEpairlEphyscLcLBosonlEphyscLcLLeptongRcOphyscLcLLeptongRsPgR);
      instance.AdoptCollectionProxyInfo(TCollectionProxyInfo::Generate(TCollectionProxyInfo::Pushback< vector<pair<phys::Boson<phys::Lepton>,phys::Lepton> > >()));
      return &instance;
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal((const vector<pair<phys::Boson<phys::Lepton>,phys::Lepton> >*)0x0); R__UseDummy(_R__UNIQUE_DICT_(Init));

   // Dictionary for non-ClassDef classes
   static TClass *vectorlEpairlEphyscLcLBosonlEphyscLcLLeptongRcOphyscLcLLeptongRsPgR_Dictionary() {
      TClass* theClass =::ROOT::GenerateInitInstanceLocal((const vector<pair<phys::Boson<phys::Lepton>,phys::Lepton> >*)0x0)->GetClass();
      vectorlEpairlEphyscLcLBosonlEphyscLcLLeptongRcOphyscLcLLeptongRsPgR_TClassManip(theClass);
   return theClass;
   }

   static void vectorlEpairlEphyscLcLBosonlEphyscLcLLeptongRcOphyscLcLLeptongRsPgR_TClassManip(TClass* ){
   }

} // end of namespace ROOT

namespace ROOT {
   // Wrappers around operator new
   static void *new_vectorlEpairlEphyscLcLBosonlEphyscLcLLeptongRcOphyscLcLLeptongRsPgR(void *p) {
      return  p ? ::new((::ROOT::Internal::TOperatorNewHelper*)p) vector<pair<phys::Boson<phys::Lepton>,phys::Lepton> > : new vector<pair<phys::Boson<phys::Lepton>,phys::Lepton> >;
   }
   static void *newArray_vectorlEpairlEphyscLcLBosonlEphyscLcLLeptongRcOphyscLcLLeptongRsPgR(Long_t nElements, void *p) {
      return p ? ::new((::ROOT::Internal::TOperatorNewHelper*)p) vector<pair<phys::Boson<phys::Lepton>,phys::Lepton> >[nElements] : new vector<pair<phys::Boson<phys::Lepton>,phys::Lepton> >[nElements];
   }
   // Wrapper around operator delete
   static void delete_vectorlEpairlEphyscLcLBosonlEphyscLcLLeptongRcOphyscLcLLeptongRsPgR(void *p) {
      delete ((vector<pair<phys::Boson<phys::Lepton>,phys::Lepton> >*)p);
   }
   static void deleteArray_vectorlEpairlEphyscLcLBosonlEphyscLcLLeptongRcOphyscLcLLeptongRsPgR(void *p) {
      delete [] ((vector<pair<phys::Boson<phys::Lepton>,phys::Lepton> >*)p);
   }
   static void destruct_vectorlEpairlEphyscLcLBosonlEphyscLcLLeptongRcOphyscLcLLeptongRsPgR(void *p) {
      typedef vector<pair<phys::Boson<phys::Lepton>,phys::Lepton> > current_t;
      ((current_t*)p)->~current_t();
   }
} // end of namespace ROOT for class vector<pair<phys::Boson<phys::Lepton>,phys::Lepton> >

namespace ROOT {
   static TClass *bitsetlE15gR_Dictionary();
   static void bitsetlE15gR_TClassManip(TClass*);
   static void *new_bitsetlE15gR(void *p = 0);
   static void *newArray_bitsetlE15gR(Long_t size, void *p);
   static void delete_bitsetlE15gR(void *p);
   static void deleteArray_bitsetlE15gR(void *p);
   static void destruct_bitsetlE15gR(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const bitset<15>*)
   {
      bitset<15> *ptr = 0;
      static ::TVirtualIsAProxy* isa_proxy = new ::TIsAProxy(typeid(bitset<15>));
      static ::ROOT::TGenericClassInfo 
         instance("bitset<15>", 2, "bitset", 747,
                  typeid(bitset<15>), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &bitsetlE15gR_Dictionary, isa_proxy, 0,
                  sizeof(bitset<15>) );
      instance.SetNew(&new_bitsetlE15gR);
      instance.SetNewArray(&newArray_bitsetlE15gR);
      instance.SetDelete(&delete_bitsetlE15gR);
      instance.SetDeleteArray(&deleteArray_bitsetlE15gR);
      instance.SetDestructor(&destruct_bitsetlE15gR);
      instance.AdoptCollectionProxyInfo(TCollectionProxyInfo::Generate(TCollectionProxyInfo::Pushback<Internal::TStdBitsetHelper< bitset<15> > >()));
      return &instance;
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal((const bitset<15>*)0x0); R__UseDummy(_R__UNIQUE_DICT_(Init));

   // Dictionary for non-ClassDef classes
   static TClass *bitsetlE15gR_Dictionary() {
      TClass* theClass =::ROOT::GenerateInitInstanceLocal((const bitset<15>*)0x0)->GetClass();
      bitsetlE15gR_TClassManip(theClass);
   return theClass;
   }

   static void bitsetlE15gR_TClassManip(TClass* ){
   }

} // end of namespace ROOT

namespace ROOT {
   // Wrappers around operator new
   static void *new_bitsetlE15gR(void *p) {
      return  p ? ::new((::ROOT::Internal::TOperatorNewHelper*)p) bitset<15> : new bitset<15>;
   }
   static void *newArray_bitsetlE15gR(Long_t nElements, void *p) {
      return p ? ::new((::ROOT::Internal::TOperatorNewHelper*)p) bitset<15>[nElements] : new bitset<15>[nElements];
   }
   // Wrapper around operator delete
   static void delete_bitsetlE15gR(void *p) {
      delete ((bitset<15>*)p);
   }
   static void deleteArray_bitsetlE15gR(void *p) {
      delete [] ((bitset<15>*)p);
   }
   static void destruct_bitsetlE15gR(void *p) {
      typedef bitset<15> current_t;
      ((current_t*)p)->~current_t();
   }
} // end of namespace ROOT for class bitset<15>

namespace {
  void TriggerDictionaryInitialization_libTreeAnalysis_Impl() {
    static const char* headers[] = {
"../..//VVXAnalysis/DataFormats/interface/DiBoson.h",
"../..//VVXAnalysis/DataFormats/interface/TypeDefs.h",
"../..//VVXAnalysis/DataFormats/interface/GenEventWeights.h",
"../..//VVXAnalysis/DataFormats/interface/MELA.h",
"../..//VVXAnalysis/DataFormats/interface/Jet.h",
"../..//VVXAnalysis/DataFormats/interface/Boson.h",
"../..//VVXAnalysis/DataFormats/interface/GenStatusBit.h",
"../..//VVXAnalysis/DataFormats/interface/Lepton.h",
"../..//VVXAnalysis/DataFormats/interface/Electron.h",
"../..//VVXAnalysis/DataFormats/interface/Particle.h",
0
    };
    static const char* includePaths[] = {
"/usr/local/root/include",
"/home/cristiano/VVXAnalysis/TreeAnalysis",
"/home/cristiano/VVXAnalysis/TreeAnalysis/../..",
"/home/cristiano/VVXAnalysis/TreeAnalysis/interface",
"/opt/local/include",
"/home/cristiano/VVXAnalysis/TreeAnalysis/../DataFormats/interface",
"/usr/local/root/include",
"/home/cristiano/VVXAnalysis/TreeAnalysis/",
0
    };
    static const char* fwdDeclCode = R"DICTFWDDCLS(
#line 1 "libTreeAnalysis dictionary forward declarations' payload"
#pragma clang diagnostic ignored "-Wkeyword-compat"
#pragma clang diagnostic ignored "-Wignored-attributes"
#pragma clang diagnostic ignored "-Wreturn-type-c-linkage"
extern int __Cling_Autoloading_Map;
namespace phys{class __attribute__((annotate("$clingAutoload$VVXAnalysis/DataFormats/interface/Lepton.h")))  Lepton;}
namespace phys{template <typename P> class __attribute__((annotate("$clingAutoload$VVXAnalysis/DataFormats/interface/Boson.h")))  Boson;
}
namespace std{template <class _T1, class _T2> struct __attribute__((annotate("$clingAutoload$bits/stl_pair.h")))  __attribute__((annotate("$clingAutoload$string")))  pair;
}
namespace std{template <typename _Tp> class __attribute__((annotate("$clingAutoload$bits/allocator.h")))  __attribute__((annotate("$clingAutoload$string")))  allocator;
}
namespace phys{class __attribute__((annotate("$clingAutoload$VVXAnalysis/DataFormats/interface/Electron.h")))  Electron;}
namespace phys{template <typename P1, typename P2> class __attribute__((annotate("$clingAutoload$VVXAnalysis/DataFormats/interface/DiBoson.h")))  DiBoson;
}
namespace phys{class __attribute__((annotate("$clingAutoload$VVXAnalysis/DataFormats/interface/Particle.h")))  Particle;}
namespace phys{class __attribute__((annotate("$clingAutoload$Jet.h")))  __attribute__((annotate("$clingAutoload$VVXAnalysis/DataFormats/interface/Lepton.h")))  Jet;}
namespace phys{class __attribute__((annotate("$clingAutoload$VVXAnalysis/DataFormats/interface/GenEventWeights.h")))  GenEventWeights;}
namespace phys{class __attribute__((annotate("$clingAutoload$VVXAnalysis/DataFormats/interface/MELA.h")))  MELA;}
)DICTFWDDCLS";
    static const char* payloadCode = R"DICTPAYLOAD(
#line 1 "libTreeAnalysis dictionary payload"

#ifndef G__VECTOR_HAS_CLASS_ITERATOR
  #define G__VECTOR_HAS_CLASS_ITERATOR 1
#endif

#define _BACKWARD_BACKWARD_WARNING_H
#include "../..//VVXAnalysis/DataFormats/interface/DiBoson.h"
#include "../..//VVXAnalysis/DataFormats/interface/TypeDefs.h"
#include "../..//VVXAnalysis/DataFormats/interface/GenEventWeights.h"
#include "../..//VVXAnalysis/DataFormats/interface/MELA.h"
#include "../..//VVXAnalysis/DataFormats/interface/Jet.h"
#include "../..//VVXAnalysis/DataFormats/interface/Boson.h"
#include "../..//VVXAnalysis/DataFormats/interface/GenStatusBit.h"
#include "../..//VVXAnalysis/DataFormats/interface/Lepton.h"
#include "../..//VVXAnalysis/DataFormats/interface/Electron.h"
#include "../..//VVXAnalysis/DataFormats/interface/Particle.h"
#include "VVXAnalysis/DataFormats/interface/Particle.h"
#include "VVXAnalysis/DataFormats/interface/Lepton.h"
#include "VVXAnalysis/DataFormats/interface/Jet.h"
#include "VVXAnalysis/DataFormats/interface/Electron.h"
#include "VVXAnalysis/DataFormats/interface/Boson.h"
#include "VVXAnalysis/DataFormats/interface/DiBoson.h"
#include "VVXAnalysis/DataFormats/interface/GenEventWeights.h"
#include "VVXAnalysis/DataFormats/interface/MELA.h"


#ifdef __CINT__

#pragma link C++ class  phys::GenEventWeights+;
#pragma link C++ class  phys::MELA+;

#pragma link C++ class  phys::Particle+;
#pragma link C++ class  phys::Lepton+;
#pragma link C++ class  phys::Jet+;
#pragma link C++ class  phys::Electron+;
#pragma link C++ class  phys::Boson<phys::Particle>+;
#pragma link C++ class  phys::Boson<phys::Lepton>+;
#pragma link C++ class  phys::Boson<phys::Electron>+;
#pragma link C++ class  phys::Boson<phys::Jet>+;
#pragma link C++ class  phys::DiBoson<phys::Particle, phys::Particle >+;
#pragma link C++ class  phys::DiBoson<phys::Lepton  , phys::Lepton >+;
#pragma link C++ class  phys::DiBoson<phys::Electron, phys::Lepton >+;
#pragma link C++ class  phys::DiBoson<phys::Lepton  , phys::Electron >+;
#pragma link C++ class  phys::DiBoson<phys::Electron, phys::Electron >+;

#pragma link C++ class  std::vector<phys::Particle>;
#pragma link C++ class  std::vector<phys::Lepton>;
#pragma link C++ class  std::vector<phys::Jet>;
#pragma link C++ class  std::vector<phys::Electron>;
#pragma link C++ class  std::vector<phys::Boson<phys::Particle> >;
#pragma link C++ class  std::vector<phys::Boson<phys::Lepton> >;
#pragma link C++ class  std::vector<phys::Boson<phys::Electron> >;
#pragma link C++ class  std::vector<phys::Boson<phys::Jet> >;
#pragma link C++ class  std::vector<phys::DiBoson<phys::Particle, phys::Particle > >;
#pragma link C++ class  std::vector<phys::DiBoson<phys::Lepton  , phys::Lepton > >;
#pragma link C++ class  std::vector<phys::DiBoson<phys::Electron, phys::Lepton > >;
#pragma link C++ class  std::vector<phys::DiBoson<phys::Lepton  , phys::Electron > >;
#pragma link C++ class  std::vector<phys::DiBoson<phys::Electron, phys::Electron > >;
#pragma link C++ class  std::pair<phys::Boson<phys::Lepton>, phys::Lepton>+;
#pragma link C++ class  std::vector<std::pair<phys::Boson<phys::Lepton>, phys::Lepton> >;


#endif

#undef  _BACKWARD_BACKWARD_WARNING_H
)DICTPAYLOAD";
    static const char* classesHeaders[]={
"phys::Boson<phys::Electron>", payloadCode, "@",
"phys::Boson<phys::Jet>", payloadCode, "@",
"phys::Boson<phys::Lepton>", payloadCode, "@",
"phys::Boson<phys::Particle>", payloadCode, "@",
"phys::DiBoson<phys::Electron,phys::Electron>", payloadCode, "@",
"phys::DiBoson<phys::Electron,phys::Lepton>", payloadCode, "@",
"phys::DiBoson<phys::Lepton,phys::Electron>", payloadCode, "@",
"phys::DiBoson<phys::Lepton,phys::Lepton>", payloadCode, "@",
"phys::DiBoson<phys::Particle,phys::Particle>", payloadCode, "@",
"phys::Electron", payloadCode, "@",
"phys::GenEventWeights", payloadCode, "@",
"phys::Jet", payloadCode, "@",
"phys::Lepton", payloadCode, "@",
"phys::MELA", payloadCode, "@",
"phys::Particle", payloadCode, "@",
nullptr};

    static bool isInitialized = false;
    if (!isInitialized) {
      TROOT::RegisterModule("libTreeAnalysis",
        headers, includePaths, payloadCode, fwdDeclCode,
        TriggerDictionaryInitialization_libTreeAnalysis_Impl, {}, classesHeaders, /*has no C++ module*/false);
      isInitialized = true;
    }
  }
  static struct DictInit {
    DictInit() {
      TriggerDictionaryInitialization_libTreeAnalysis_Impl();
    }
  } __TheDictionaryInitializer;
}
void TriggerDictionaryInitialization_libTreeAnalysis() {
  TriggerDictionaryInitialization_libTreeAnalysis_Impl();
}
