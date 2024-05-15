#ifndef CHANNELDATA_CRPTCR_CD2_HH
#define CHANNELDATA_CRPTCR_CD2_HH

#include "TObject.h"

class ChannelData : public TObject
{
  public:
    UShort_t fModule;
    UShort_t fChannel;
    UShort_t fEnergy;
    Long64_t fTS;
    UShort_t fTSGroup;
    Long64_t fTSFull;

    ChannelData() { Clear(); }
    ~ChannelData() {}

    void Clear(Option_t* option="") {
      fModule = -1;
      fChannel = -1;
      fEnergy = -1;
      fTS = -1;
      fTSGroup = -1;
      fTSFull = -1;
    }

    void SetData(UShort_t module, UShort_t channel, UShort_t energy,
                 Long64_t ts, UShort_t tsGroup, Long64_t tsFull)
    {
      fModule = module;
      fChannel = channel;;
      fEnergy = energy;;
      fTS = ts;
      fTSGroup = tsGroup;
      fTSFull = tsFull;
    }

  ClassDef(ChannelData,1)
};

#endif
