#ifndef CHANNELDATA_CRPTCR_CD2_HH
#define CHANNELDATA_CRPTCR_CD2_HH

#include "TObject.h"

class ChannelData : public TObject
{
  public:
    UShort_t fModule;
    UShort_t fChannel;
    UShort_t fEnergy;
    Long64_t fTimeStamp0;
    UShort_t fTSGroup;
    Long64_t fTimeStamp;

    ChannelData() { Clear(); }
    ~ChannelData() {}

    void Clear(Option_t* option="") {
      fModule = -1;
      fChannel = -1;
      fEnergy = -1;
      fTimeStamp0 = -1;
      fTSGroup = -1;
      fTimeStamp = -1;
    }

    void SetData(UShort_t module, UShort_t channel, UShort_t energy,
                 Long64_t ts0, UShort_t tsGroup, Long64_t ts)
    {
      fModule = module;
      fChannel = channel;
      fEnergy = energy;
      fTimeStamp0 = ts0;
      fTSGroup = tsGroup;
      fTimeStamp = ts;
    }

  ClassDef(ChannelData,1)
};

#endif
