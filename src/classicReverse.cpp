#define STRICT_R_HEADER
#include <RcppArmadillo.h>

// EVID = 0; Observations
// EVID = 1; is illegal, but converted from NONMEM
// EVID = 2; Non-observation, possibly covariate
// EVID = 3; Reset ODE states to zero; Non-observation event
// EVID = 4; Reset and then dose event;  Illegal
// EVID = 9; Non-observation event to ini system at time zero; This is to set the INIs at the correct place.
// EVID = 10-99; mtime events (from ODE system)
// When EVID > 100
// EVID: ## # ## ##
//       c2 I c1 xx
// c2 = Compartment numbers over 100
//  I = Infusion Flag/ Special event flag
#define EVIDF_NORMAL 0

#define EVIDF_INF_RATE 1
#define EVIDF_INF_DUR  2

#define EVIDF_REPLACE  4
#define EVIDF_MULT     5

#define EVIDF_MODEL_DUR_ON   8
#define EVIDF_MODEL_DUR_OFF  6

#define EVIDF_MODEL_RATE_ON  9
#define EVIDF_MODEL_RATE_OFF 7
//      0 = no Infusion
//      1 = Infusion, AMT=rate (mg/hr for instance)
//      2 = Infusion, duration is fixed
//      4 = Replacement event
//      5 = Multiplication event
//      6 = Turn off modeled duration
//      7 = Turn off modeled rate compartment
//      8 = Duration is modeled, AMT=dose; Rate = AMT/(Modeled Duration) NONMEM RATE=-2
//      9 = Rate is modeled, AMT=dose; Duration = AMT/(Modeled Rate) NONMEM RATE=-1
// c1 = Compartment numbers below 99
// xx = 1, regular event
// xx = 10, steady state event SS=1
// xx = 20, steady state event + last observed info.
// xx = 30, Turn off compartment
// xx = 40, Steady state constant infusion
// xx = 50, Phantom event, used for transit compartments
// Steady state events need a II data item > 0
#define EVID0_REGULAR 1
#define EVID0_SS 10
#define EVID0_SS2 20
#define EVID0_OFF 30
#define EVID0_SSINF 40 
#define EVID0_PHANTOM 50

static inline void getWh(int evid, int *wh, int *cmt, int *wh100, int *whI, int *wh0,
                         int linNcmt, int linKa, int neq){
  *wh = evid;
  *cmt = 0;
  *wh100 = std::floor(*wh/1e5L);
  *whI   = std::floor(*wh/1e4L-*wh100*10);
  *wh    = *wh - *wh100*1e5 - (*whI-1)*1e4;
  *wh0 = std::floor((*wh%10000)/100);
  *cmt = *wh0 - 1 + *wh100*100;
  *wh0 = evid - *wh100*1e5 - *whI*1e4 - *wh0*100;
  if (linNcmt != 0) {
    if (linKa) {
      switch (*cmt) {
      case 0:
				*cmt = neq;
				break;
      case 1:
				*cmt = neq+1;
				break;
      case 2:
				*cmt -= 2;
				break;
      }
    } else {
      if (*cmt == 0) {
				*cmt = neq;
      } else {
				*cmt -= 1;
      }
    }
  }
}

static inline int getSs(int wh0) {
  if (wh0 == EVID0_SS2) {
    return  2;
  } else if (wh0 == EVID0_SS) {
    return 1;
  } else if (wh0 == EVID0_SSINF) {
    // corresponds to ss=1, amt=0, ii=0
    // both amt and ii should be 0, and rate should contain the infusion rate
    return 1;
  }
  return 0;
}

using namespace Rcpp;

//[[Rcpp::export]]
List convertDataBack(IntegerVector id, NumericVector time, NumericVector amt, NumericVector ii,
                     IntegerVector evid, IntegerVector cmt,
                     int cmtOffset=0, int linNcmt=0, int linKa=0, int neq=0) {
  IntegerVector newEvid(evid.size());
  IntegerVector newSs(evid.size());
  IntegerVector newDvid(evid.size());
  LogicalVector keepItem(evid.size());
  LogicalVector newAmt(evid.size());
  NumericVector newRate(evid.size());
  NumericVector newTinf(evid.size());
  IntegerVector newCmt(evid.size());
  int wh=0, cmt0=0, wh100=0, whI=0, wh0=0;
  bool turnOffCmt=false;
  bool hasTinf = false;
  bool hasRate = false;
  bool hasPhantom = false;
  double curAmt=0.0;
  for (unsigned int i = 0; i < evid.size(); ++i) {
    int curEvid = evid[i];
    // put in defaults
    newEvid[i] = evid[i];
    newSs[i] = 0;
    newDvid[i] = 0;
    newRate[i] = 0;
    newTinf[i] = 0;
    newAmt[i] = amt[i];
    newCmt[i] = cmt[i];
    if (curEvid == 0 || curEvid==2 || curEvid == 3) {
      // 0, 2 and 3 are preserved
      newDvid[i] = curEvid == 3 ? 0 : cmt[i]-cmtOffset ;
      keepItem[i] = true;
    } else if ((curEvid) >= 9 && (curEvid) <= 99) {
      // mtimes and zero events are ignored
      keepItem[i] = false;
    } else {
      // these are doses
      getWh(curEvid, &wh, &cmt0, &wh100, &whI, &wh0,
            linNcmt, linKa, neq);
      if (wh0 ==EVID0_OFF) {
        // turn off a compartment; supported in nonmem
        newCmt[i] = -cmt[i];
        turnOffCmt=true;
        keepItem[i] = true;
      } else if (wh0 == EVID0_PHANTOM) {
        keepItem[i] = false;
        hasPhantom = true;
      } else {
        switch (whI) {
        case EVIDF_MODEL_RATE_ON: // modeled rate.
          newEvid[i] = 1;
          newSs[i] =getSs(wh0);
          newRate[i] = -1;
          newAmt[i] = amt[i];
          keepItem[i] = true;
          break;
        case EVIDF_MODEL_DUR_ON: // modeled duration.
          newEvid[i] = 1;
          newSs[i] = getSs(wh0);
          newRate[i] = -2;
          newAmt[i] = amt[i];
          keepItem[i] = true;
          break;
        case EVIDF_MODEL_RATE_OFF: // End modeled rate
        case EVIDF_MODEL_DUR_OFF: // end modeled duration
          keepItem[i] = false;
          break;
        case EVIDF_INF_DUR:
        case EVIDF_INF_RATE:
          // classic rxode2 uses rate instead of amount here
          curAmt = amt[i];
          if (amt[i] > 0) {
            bool foundOff=false;
            double t2=0;
            int curId = id[i];
            for (int j=i; j < evid.size(); j++){
              if (id[j] != curId) {
                break;
              }
              if (evid[j] == curEvid && curAmt == -amt[j]) {
                t2 = time[j];
                foundOff = true;
                break;
              }
            }
            if (foundOff) {
              double dur = t2-time[i];
              double rate = amt[i];
              curAmt = rate*dur;
              if (whI==EVIDF_INF_DUR) {
                newTinf[i] = dur;
                hasTinf = true;
              } else {
                newRate[i] = rate;
                hasRate = true;
              }
              newEvid[i] = 1;
              newSs[i] = getSs(wh0);
              newAmt[i] = curAmt;
              keepItem[i] = true;
            } else {
              keepItem[i] = false;
            }
          } else {
            keepItem[i] = false;
          }
        case EVIDF_REPLACE:
          keepItem[i] = false;
        case EVIDF_MULT:
          keepItem[i] = false;
        case EVIDF_NORMAL:
          newEvid[i] = 1;
          newSs[i] =getSs(wh0);
          newAmt[i] = amt[i];
          keepItem[i] = true;          
        }
      }
    }
  }
  // now return the dataset
  return List::create(_["df"]=DataFrame::create(
                                                _["EVID"]=newEvid,
                                                _["SS"]=newSs,
                                                _["DVID"]=newDvid,
                                                _["AMT"]=newAmt,
                                                _["RATE"]=newRate,
                                                _["TINF"]=newTinf,
                                                _["CMT"]=newCmt,
                                                _[".nlmixrKeep"]=keepItem
                                                ),
                      _["turnOffCmt"]=turnOffCmt,
                      _["hasTinf"]=hasTinf,
                      _["hasRate"]=hasRate,
                      _["hasPhantom"]=hasPhantom);
}
