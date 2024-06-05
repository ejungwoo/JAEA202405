#ifndef DETECTOR_SETTING_HH
#define DETECTOR_SETTING_HH

#ifndef NUMBER_OF_MODULES
#define NUMBER_OF_MODULES 7
#define NUMBER_OF_CHANNELS 16
#define NUMBER_OF_DETECTORS 7
#define NUM_MAX_DCH 32
const int kXX  = 0;
const int kdE  = 1;
const int kS1J = 2;
const int kS3J = 3;
const int kS3O = 4;
const int kSC  = 5;
const int kFC  = 6;
#endif

#include "TString.h"
#include <fstream>
#include <iostream>
using namespace std;

#include "TObject.h"

class DetectorSetting : public TObject
{
    public:
        DetectorSetting();
        ~DetectorSetting() {}

        int  FENumber_MIdx   (int FENumber) { return fMapFEToModuleIndex[FENumber]; }

        int  MIdxMCh_GID     (int midx, int mch) { return (midx*fNumChannels + mch); }
        int  MIdxMCh_Det     (int midx, int mch) { return fMapDetectorType[midx][mch]; }
        int  MIdxMCh_DCh     (int midx, int mch) { return fMapDetectorChannel[midx][mch]; }
        int  MIdxMCh_Replaced(int midx, int mch) { return fMapDetectorReplaced[midx][mch]; }
        int  MIdxMCh_Group   (int midx, int mch) { return fMapDetectorGroup[midx][mch]; }
        int  MIdxMCh_RFMod   (int midx, int mch) { return fMapDetectorRFMod[midx][mch]; }
        int  MIdxMCh_RFMCh   (int midx, int mch) { return fMapDetectorRFMCh[midx][mch]; }

        int  GID_MIdx  (int gid) { return fMapGlobalIDToMIdx[gid]; }
        int  GID_MCh   (int gid) { return fMapGlobalIDToMCh[gid]; }
        int  GID_Det   (int gid) { return fMapGlobalIDToDet[gid]; }
        int  GID_DCh   (int gid) { return fMapGlobalIDToDCh[gid]; }

        int  DetDCh_GID  (int det, int dch) { return fMapDetectorToGlobalID[det][dch]; }

        double  S1Ch_Angle    (int dch)  { return fMapS1ChToAngle[dch];  }
        double  S1Ch_Angle1   (int dch)  { return fMapS1ChToAngle1[dch]; }
        double  S1Ch_Angle2   (int dch)  { return fMapS1ChToAngle2[dch]; }
        int     S1Ch_Strip    (int dch)  { return fMapS1ChToStrip[dch];  }
        double  S3Ch_Angle    (int dch)  { return fMapS3ChToAngle[dch];  }
        double  S3Ch_Angle1   (int dch)  { return fMapS3ChToAngle1[dch]; }
        double  S3Ch_Angle2   (int dch)  { return fMapS3ChToAngle2[dch]; }
        int     S3Ch_Strip    (int dch)  { return fMapS3ChToStrip[dch];  }

        TString GetDetectorName(int det) { return fDetectorName[det]; }

        TString GetDetectorTitle(int midx, int mch, bool addChannel) {
            int det = fMapDetectorType[midx][mch];
            int dch = fMapDetectorChannel[midx][mch];
            TString detectorTitle = fDetectorName[det];
            if (addChannel)
                detectorTitle = detectorTitle + Form("-%d",dch);
            return detectorTitle;
        }

    private:
        const int fNumModules = NUMBER_OF_MODULES;
        const int fNumChannels = NUMBER_OF_CHANNELS;
        const int fNumDetectors = NUMBER_OF_DETECTORS;
        const int fMaxDCh = NUM_MAX_DCH;
        const int fNumCh = NUMBER_OF_MODULES*NUMBER_OF_CHANNELS;
        TString fDetectorName[NUMBER_OF_DETECTORS];

        int    fMapFEToModuleIndex[20];
        int  **fMapDetectorType;
        int  **fMapDetectorChannel;
        bool **fMapDetectorReplaced;
        int  **fMapDetectorGroup;
        int  **fMapDetectorRFMod;
        int  **fMapDetectorRFMCh;
        int  **fMapDetectorToGlobalID;
        int   *fMapGlobalIDToMIdx;
        int   *fMapGlobalIDToMCh;
        int   *fMapGlobalIDToDet;
        int   *fMapGlobalIDToDCh;
        double fMapS1ChToAngle[33];
        double fMapS1ChToAngle1[33];
        double fMapS1ChToAngle2[33];
        double fMapS3ChToAngle1[33];
        double fMapS3ChToAngle[33];
        double fMapS3ChToAngle2[33];
        int    fMapS1ChToStrip[33];
        int    fMapS3ChToStrip[33];

    private:
        TString fMapFileName = "source/ModCh.in";
        TString fDetFileName = "source/MapDetector.in";

    public:
        static DetectorSetting* GetDetectorSetting();

    private:
        static DetectorSetting* fInstance;

    ClassDef(DetectorSetting,1)
};

//DetectorSetting *getDet() { return DetectorSetting::GetDetectorSetting(); }

#endif
