#ifndef CHANNELDATA_CRPTCR_CD2_HH
#define CHANNELDATA_CRPTCR_CD2_HH

#define WRITE_ALL_DATA

#include "TObject.h"

class ChannelData : public TObject
{
  public:
#ifdef WRITE_ALL_DATA
    UShort_t module;
    UShort_t id;
    Long64_t ts0;
    UShort_t tsG;
    Long64_t ts;
#endif
    UShort_t gid;
    UShort_t energy;

    ChannelData() { Clear(); }
    ~ChannelData() {}

    void Clear(Option_t* option="") {
#ifdef WRITE_ALL_DATA
      module = -1;
      id = -1;
      ts0 = -1;
      tsG = -1;
      ts = -1;
#endif
      gid = -1;
      energy = -1;
    }

    void SetData(UShort_t module_, UShort_t id_, Int_t gid_, UShort_t energy_, Long64_t ts0_, UShort_t tsG_, Long64_t ts_)
    {
#ifdef WRITE_ALL_DATA
      module = module_;
      id = id_;
      ts0 = ts0_;
      tsG = tsG_;
      ts = ts_;
#endif
      gid = gid_;
      energy = energy_;
    }

  ClassDef(ChannelData,1)
};

#endif
