// Do NOT change. Changes will be lost next time file is generated

#define R__DICTIONARY_FILENAME hit_hits_class

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
#include "hit_hits_class.cc"
#include "hit_hits_class.h"

// Header files passed via #pragma extra_include

namespace ROOT {
   static void *new_hit(void *p = 0);
   static void *newArray_hit(Long_t size, void *p);
   static void delete_hit(void *p);
   static void deleteArray_hit(void *p);
   static void destruct_hit(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const ::hit*)
   {
      ::hit *ptr = 0;
      static ::TVirtualIsAProxy* isa_proxy = new ::TInstrumentedIsAProxy< ::hit >(0);
      static ::ROOT::TGenericClassInfo 
         instance("hit", ::hit::Class_Version(), "hit_hits_class.h", 10,
                  typeid(::hit), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &::hit::Dictionary, isa_proxy, 4,
                  sizeof(::hit) );
      instance.SetNew(&new_hit);
      instance.SetNewArray(&newArray_hit);
      instance.SetDelete(&delete_hit);
      instance.SetDeleteArray(&deleteArray_hit);
      instance.SetDestructor(&destruct_hit);
      return &instance;
   }
   TGenericClassInfo *GenerateInitInstance(const ::hit*)
   {
      return GenerateInitInstanceLocal((::hit*)0);
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_(Init) = GenerateInitInstanceLocal((const ::hit*)0x0); R__UseDummy(_R__UNIQUE_(Init));
} // end of namespace ROOT

namespace ROOT {
   static void *new_hitcollection(void *p = 0);
   static void *newArray_hitcollection(Long_t size, void *p);
   static void delete_hitcollection(void *p);
   static void deleteArray_hitcollection(void *p);
   static void destruct_hitcollection(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const ::hitcollection*)
   {
      ::hitcollection *ptr = 0;
      static ::TVirtualIsAProxy* isa_proxy = new ::TInstrumentedIsAProxy< ::hitcollection >(0);
      static ::ROOT::TGenericClassInfo 
         instance("hitcollection", ::hitcollection::Class_Version(), "hit_hits_class.h", 37,
                  typeid(::hitcollection), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &::hitcollection::Dictionary, isa_proxy, 4,
                  sizeof(::hitcollection) );
      instance.SetNew(&new_hitcollection);
      instance.SetNewArray(&newArray_hitcollection);
      instance.SetDelete(&delete_hitcollection);
      instance.SetDeleteArray(&deleteArray_hitcollection);
      instance.SetDestructor(&destruct_hitcollection);
      return &instance;
   }
   TGenericClassInfo *GenerateInitInstance(const ::hitcollection*)
   {
      return GenerateInitInstanceLocal((::hitcollection*)0);
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_(Init) = GenerateInitInstanceLocal((const ::hitcollection*)0x0); R__UseDummy(_R__UNIQUE_(Init));
} // end of namespace ROOT

//______________________________________________________________________________
atomic_TClass_ptr hit::fgIsA(0);  // static to hold class pointer

//______________________________________________________________________________
const char *hit::Class_Name()
{
   return "hit";
}

//______________________________________________________________________________
const char *hit::ImplFileName()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::hit*)0x0)->GetImplFileName();
}

//______________________________________________________________________________
int hit::ImplFileLine()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::hit*)0x0)->GetImplFileLine();
}

//______________________________________________________________________________
TClass *hit::Dictionary()
{
   fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::hit*)0x0)->GetClass();
   return fgIsA;
}

//______________________________________________________________________________
TClass *hit::Class()
{
   if (!fgIsA.load()) { R__LOCKGUARD2(gInterpreterMutex); fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::hit*)0x0)->GetClass(); }
   return fgIsA;
}

//______________________________________________________________________________
atomic_TClass_ptr hitcollection::fgIsA(0);  // static to hold class pointer

//______________________________________________________________________________
const char *hitcollection::Class_Name()
{
   return "hitcollection";
}

//______________________________________________________________________________
const char *hitcollection::ImplFileName()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::hitcollection*)0x0)->GetImplFileName();
}

//______________________________________________________________________________
int hitcollection::ImplFileLine()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::hitcollection*)0x0)->GetImplFileLine();
}

//______________________________________________________________________________
TClass *hitcollection::Dictionary()
{
   fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::hitcollection*)0x0)->GetClass();
   return fgIsA;
}

//______________________________________________________________________________
TClass *hitcollection::Class()
{
   if (!fgIsA.load()) { R__LOCKGUARD2(gInterpreterMutex); fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::hitcollection*)0x0)->GetClass(); }
   return fgIsA;
}

//______________________________________________________________________________
void hit::Streamer(TBuffer &R__b)
{
   // Stream an object of class hit.

   if (R__b.IsReading()) {
      R__b.ReadClassBuffer(hit::Class(),this);
   } else {
      R__b.WriteClassBuffer(hit::Class(),this);
   }
}

namespace ROOT {
   // Wrappers around operator new
   static void *new_hit(void *p) {
      return  p ? new(p) ::hit : new ::hit;
   }
   static void *newArray_hit(Long_t nElements, void *p) {
      return p ? new(p) ::hit[nElements] : new ::hit[nElements];
   }
   // Wrapper around operator delete
   static void delete_hit(void *p) {
      delete ((::hit*)p);
   }
   static void deleteArray_hit(void *p) {
      delete [] ((::hit*)p);
   }
   static void destruct_hit(void *p) {
      typedef ::hit current_t;
      ((current_t*)p)->~current_t();
   }
} // end of namespace ROOT for class ::hit

//______________________________________________________________________________
void hitcollection::Streamer(TBuffer &R__b)
{
   // Stream an object of class hitcollection.

   if (R__b.IsReading()) {
      R__b.ReadClassBuffer(hitcollection::Class(),this);
   } else {
      R__b.WriteClassBuffer(hitcollection::Class(),this);
   }
}

namespace ROOT {
   // Wrappers around operator new
   static void *new_hitcollection(void *p) {
      return  p ? new(p) ::hitcollection : new ::hitcollection;
   }
   static void *newArray_hitcollection(Long_t nElements, void *p) {
      return p ? new(p) ::hitcollection[nElements] : new ::hitcollection[nElements];
   }
   // Wrapper around operator delete
   static void delete_hitcollection(void *p) {
      delete ((::hitcollection*)p);
   }
   static void deleteArray_hitcollection(void *p) {
      delete [] ((::hitcollection*)p);
   }
   static void destruct_hitcollection(void *p) {
      typedef ::hitcollection current_t;
      ((current_t*)p)->~current_t();
   }
} // end of namespace ROOT for class ::hitcollection

namespace ROOT {
   static TClass *vectorlEhitgR_Dictionary();
   static void vectorlEhitgR_TClassManip(TClass*);
   static void *new_vectorlEhitgR(void *p = 0);
   static void *newArray_vectorlEhitgR(Long_t size, void *p);
   static void delete_vectorlEhitgR(void *p);
   static void deleteArray_vectorlEhitgR(void *p);
   static void destruct_vectorlEhitgR(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const vector<hit>*)
   {
      vector<hit> *ptr = 0;
      static ::TVirtualIsAProxy* isa_proxy = new ::TIsAProxy(typeid(vector<hit>));
      static ::ROOT::TGenericClassInfo 
         instance("vector<hit>", -2, "vector", 214,
                  typeid(vector<hit>), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &vectorlEhitgR_Dictionary, isa_proxy, 4,
                  sizeof(vector<hit>) );
      instance.SetNew(&new_vectorlEhitgR);
      instance.SetNewArray(&newArray_vectorlEhitgR);
      instance.SetDelete(&delete_vectorlEhitgR);
      instance.SetDeleteArray(&deleteArray_vectorlEhitgR);
      instance.SetDestructor(&destruct_vectorlEhitgR);
      instance.AdoptCollectionProxyInfo(TCollectionProxyInfo::Generate(TCollectionProxyInfo::Pushback< vector<hit> >()));
      return &instance;
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_(Init) = GenerateInitInstanceLocal((const vector<hit>*)0x0); R__UseDummy(_R__UNIQUE_(Init));

   // Dictionary for non-ClassDef classes
   static TClass *vectorlEhitgR_Dictionary() {
      TClass* theClass =::ROOT::GenerateInitInstanceLocal((const vector<hit>*)0x0)->GetClass();
      vectorlEhitgR_TClassManip(theClass);
   return theClass;
   }

   static void vectorlEhitgR_TClassManip(TClass* ){
   }

} // end of namespace ROOT

namespace ROOT {
   // Wrappers around operator new
   static void *new_vectorlEhitgR(void *p) {
      return  p ? ::new((::ROOT::Internal::TOperatorNewHelper*)p) vector<hit> : new vector<hit>;
   }
   static void *newArray_vectorlEhitgR(Long_t nElements, void *p) {
      return p ? ::new((::ROOT::Internal::TOperatorNewHelper*)p) vector<hit>[nElements] : new vector<hit>[nElements];
   }
   // Wrapper around operator delete
   static void delete_vectorlEhitgR(void *p) {
      delete ((vector<hit>*)p);
   }
   static void deleteArray_vectorlEhitgR(void *p) {
      delete [] ((vector<hit>*)p);
   }
   static void destruct_vectorlEhitgR(void *p) {
      typedef vector<hit> current_t;
      ((current_t*)p)->~current_t();
   }
} // end of namespace ROOT for class vector<hit>

namespace {
  void TriggerDictionaryInitialization_hit_hits_class_Impl() {
    static const char* headers[] = {
"hit_hits_class.cc",
"hit_hits_class.h",
0
    };
    static const char* includePaths[] = {
"/usr/local/include/root",
"/home/playgame555/my_RPI_ana/Ntuplizer/",
0
    };
    static const char* fwdDeclCode = R"DICTFWDDCLS(
#line 1 "hit_hits_class dictionary forward declarations' payload"
#pragma clang diagnostic ignored "-Wkeyword-compat"
#pragma clang diagnostic ignored "-Wignored-attributes"
#pragma clang diagnostic ignored "-Wreturn-type-c-linkage"
extern int __Cling_Autoloading_Map;
class __attribute__((annotate("$clingAutoload$hit_hits_class.h")))  __attribute__((annotate("$clingAutoload$hit_hits_class.cc")))  hit;
namespace std{template <typename _Tp> class __attribute__((annotate("$clingAutoload$bits/allocator.h")))  __attribute__((annotate("$clingAutoload$string")))  allocator;
}
class __attribute__((annotate("$clingAutoload$hit_hits_class.h")))  __attribute__((annotate("$clingAutoload$hit_hits_class.cc")))  hitcollection;
)DICTFWDDCLS";
    static const char* payloadCode = R"DICTPAYLOAD(
#line 1 "hit_hits_class dictionary payload"

#ifndef G__VECTOR_HAS_CLASS_ITERATOR
  #define G__VECTOR_HAS_CLASS_ITERATOR 1
#endif

#define _BACKWARD_BACKWARD_WARNING_H
#include "hit_hits_class.cc"
#include "hit_hits_class.h"

#undef  _BACKWARD_BACKWARD_WARNING_H
)DICTPAYLOAD";
    static const char* classesHeaders[]={
"hit", payloadCode, "@",
"hitcollection", payloadCode, "@",
nullptr};

    static bool isInitialized = false;
    if (!isInitialized) {
      TROOT::RegisterModule("hit_hits_class",
        headers, includePaths, payloadCode, fwdDeclCode,
        TriggerDictionaryInitialization_hit_hits_class_Impl, {}, classesHeaders);
      isInitialized = true;
    }
  }
  static struct DictInit {
    DictInit() {
      TriggerDictionaryInitialization_hit_hits_class_Impl();
    }
  } __TheDictionaryInitializer;
}
void TriggerDictionaryInitialization_hit_hits_class() {
  TriggerDictionaryInitialization_hit_hits_class_Impl();
}
