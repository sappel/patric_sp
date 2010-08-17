#ifndef _BUMP_H_
#define _BUMP_H_

class SectorMap;
class BeamLine;

class Bump{
 private:
  double decr[4];
  SectorMap *kickers[4];

 public:
  Bump(BeamLine* bl, double emit_x, double rmsToFull, double dp0, double Q_x,
       double &offset, double x_septum, double inj_angle, unsigned &max_inj,
       double inj_phase, unsigned sept_return, int id);
  void decrement();
};

#endif  // _BUMP_H_
