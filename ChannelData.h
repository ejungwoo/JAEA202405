#ifndef CHANNELDATA_CRPTCR_CD2_HH
#define CHANNELDATA_CRPTCR_CD2_HH

#define USE_RICH_DATA

#include "TObject.h"

class ChannelData : public TObject
{
  public:
#ifdef USE_RICH_DATA
    Short_t module;
    Short_t id;
    Long64_t ts0;
    Short_t tsG;
    Long64_t ts;
#endif
    Short_t gid;
    Short_t adc;
    Double_t energy;

    ChannelData() { Clear(); }
    ~ChannelData() {}

    void Clear(Option_t* option="") {
#ifdef USE_RICH_DATA
      module = 0;
      id = 0;
      ts0 = 0;
      tsG = 0;
      ts = 0;
#endif
      gid = 0;
      adc = 0;
      energy = 0.;
    }

    void SetData(Short_t module_, Short_t id_, Int_t gid_, Short_t adc_, Double_t energy_, Long64_t ts0_, Short_t tsG_, Long64_t ts_)
    {
#ifdef USE_RICH_DATA
      module = module_;
      id = id_;
      ts0 = ts0_;
      tsG = tsG_;
      ts = ts_;
#endif
      gid = gid_;
      adc = adc_;
      energy = energy_;
    }

  ClassDef(ChannelData,1)
};

#endif
