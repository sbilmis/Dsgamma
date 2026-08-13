<|"TraceAA" -> (-2*(mQ^2 + 2*mQ*ms + ms^2 - s)*(mQ^2 - 2*mQ*ms + ms^2 + 2*s))/
   (3*s), "TraceAB" -> (2*(mQ - ms)*(-mQ - ms + Sqrt[s])*(mQ + ms + Sqrt[s])*
    Conjugate[Sqrt[s]])/(Sqrt[s]*Conjugate[mQ + ms]), 
 "TraceBB" -> (2*(-mQ - ms + Sqrt[s])*(mQ + ms + Sqrt[s])*
    (2*mQ^2 - 4*mQ*ms + 2*ms^2 + s)*Conjugate[Sqrt[s]])/
   (3*Sqrt[s]*Abs[mQ + ms]^2), 
 "rhoAA" -> ((-(mQ + ms)^2 + s)*Sqrt[(-(mQ - ms)^2 + s)*(-(mQ + ms)^2 + s)]*
    ((mQ - ms)^2 + 2*s))/(8*Pi^2*s^2), 
 "rhoAB" -> (3*(mQ - ms)*(-(mQ + ms)^2 + s)*
    Sqrt[(-(mQ - ms)^2 + s)*(-(mQ + ms)^2 + s)])/(8*(mQ + ms)*Pi^2*s), 
 "rhoBA" -> (3*(mQ - ms)*(-(mQ + ms)^2 + s)*
    Sqrt[(-(mQ - ms)^2 + s)*(-(mQ + ms)^2 + s)])/(8*(mQ + ms)*Pi^2*s), 
 "rhoBB" -> ((2*(mQ - ms)^2 + s)*(-(mQ + ms)^2 + s)*
    Sqrt[(-(mQ - ms)^2 + s)*(-(mQ + ms)^2 + s)])/(8*(mQ + ms)^2*Pi^2*s), 
 "SpectralDeterminant" -> ((mQ^2 - 2*mQ*ms + ms^2 - s)^3*
    (mQ^2 + 2*mQ*ms + ms^2 - s)^3)/(32*(mQ + ms)^2*Pi^4*s^3), 
 "d3Matrix" -> {{(mQ*ss)/E^(mQ^2/M2), (mQ^2*ss)/(E^(mQ^2/M2)*(mQ + ms))}, 
   {(mQ^2*ss)/(E^(mQ^2/M2)*(mQ + ms)), (mQ^3*ss)/(E^(mQ^2/M2)*(mQ + ms)^2)}}, 
 "d3msMatrix" -> {{((-1/2*(mQ^2*ms)/M2 + (mQ^3*ms^2)/(2*M2^2))*ss)/
     E^(mQ^2/M2), (mQ*(-1/2*(mQ^2*ms)/M2 + (mQ^3*ms^2)/(2*M2^2))*ss)/
     (E^(mQ^2/M2)*(mQ + ms))}, 
   {(mQ*(-1/2*(mQ^2*ms)/M2 + (mQ^3*ms^2)/(2*M2^2))*ss)/
     (E^(mQ^2/M2)*(mQ + ms)), 
    (mQ^2*(-1/2*(mQ^2*ms)/M2 + (mQ^3*ms^2)/(2*M2^2))*ss)/
     (E^(mQ^2/M2)*(mQ + ms)^2)}}, "d5Matrix" -> 
  {{-1/4*(mixedSS*mQ^3)/(E^(mQ^2/M2)*M2^2), 
    -1/4*(mixedSS*(-1 - mQ^2/M2 + mQ^4/M2^2))/(E^(mQ^2/M2)*(mQ + ms))}, 
   {-1/4*(mixedSS*(-1 - mQ^2/M2 + mQ^4/M2^2))/(E^(mQ^2/M2)*(mQ + ms)), 
    (mixedSS*(mQ^3/(2*M2) - mQ^5/(4*M2^2)))/(E^(mQ^2/M2)*(mQ + ms)^2)}}, 
 "PhysicalMatrix" -> {{piBB*Cos[theta]^2 + 2*piAB*Cos[theta]*Sin[theta] + 
     piAA*Sin[theta]^2, piAB*Cos[theta]^2 + piAA*Cos[theta]*Sin[theta] - 
     piBB*Cos[theta]*Sin[theta] - piAB*Sin[theta]^2}, 
   {piAB*Cos[theta]^2 + piAA*Cos[theta]*Sin[theta] - 
     piBB*Cos[theta]*Sin[theta] - piAB*Sin[theta]^2, 
    piAA*Cos[theta]^2 - 2*piAB*Cos[theta]*Sin[theta] + piBB*Sin[theta]^2}}, 
 "Twist2PreBorelRatio" -> k2/(mQ*(mQ + ms)), "Twist2PostBorelRatio" -> 
  mQ/(mQ + ms)|>
