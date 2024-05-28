#ifndef CHANNELDATA_CRPTCR_CD2_HH
#define CHANNELDATA_CRPTCR_CD2_HH

//#define USE_RICH_DATA
#define USE_DETECTOR_IDS

#include "TObject.h"

class ChannelData : public TObject
{
    public:
#ifdef USE_RICH_DATA
        Long64_t ts0;
        Short_t tsG;
#endif
#ifdef USE_DETECTOR_IDS
        Short_t det;
        Short_t dch;
#endif
        Short_t midx;
        Short_t mch;
        Long64_t ts;
        Short_t gid;
        Short_t adc;
        Double_t energy;

        ChannelData() { Clear(); }
        ~ChannelData() {}

        void Clear(Option_t* option="") {
#ifdef USE_RICH_DATA
            ts0 = -1;
            tsG = -1;
#endif
#ifdef USE_DETECTOR_IDS
            det = -1;
            dch = -1;
#endif
            midx = -1;
            mch = -1;
            ts = -1;
            gid = -1;
            adc = -1;
            energy = -1.;
        }

        void SetData(Short_t midx_, Short_t mch_, Short_t det_, Short_t dch_, Int_t gid_, Short_t adc_, Double_t energy_, Long64_t ts0_, Short_t tsG_, Long64_t ts_)
        {
#ifdef USE_RICH_DATA
            ts0 = ts0_;
            tsG = tsG_;
#endif
#ifdef USE_DETECTOR_IDS
            det = det_;
            dch = dch_;
#endif
            midx = midx_;
            mch = mch_;
            ts = ts_;
            gid = gid_;
            adc = adc_;
            energy = energy_;
        }

        ClassDef(ChannelData,1)
};

#endif
