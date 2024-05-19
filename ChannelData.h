#ifndef CHANNELDATA_CRPTCR_CD2_HH
#define CHANNELDATA_CRPTCR_CD2_HH

//#define USE_RICH_DATA

#include "TObject.h"

class ChannelData : public TObject
{
  public:
#ifdef USE_RICH_DATA
    Long64_t ts0;
    Short_t tsG;
#endif
    Short_t midx;
    Short_t id;
    Long64_t ts;
    Short_t gid;
    Short_t adc;
    Double_t energy;

    ChannelData() { Clear(); }
    ~ChannelData() {}

    void Clear(Option_t* option="") {
#ifdef USE_RICH_DATA
      ts0 = 0;
      tsG = 0;
#endif
      midx = 0;
      id = 0;
      ts = 0;
      gid = 0;
      adc = 0;
      energy = 0.;
    }

    void SetData(Short_t midx_, Short_t id_, Int_t gid_, Short_t adc_, Double_t energy_, Long64_t ts0_, Short_t tsG_, Long64_t ts_)
    {
#ifdef USE_RICH_DATA
      ts0 = ts0_;
      tsG = tsG_;
#endif
      midx = midx_;
      id = id_;
      ts = ts_;
      gid = gid_;
      adc = adc_;
      energy = energy_;
    }

  ClassDef(ChannelData,1)
};

#endif
