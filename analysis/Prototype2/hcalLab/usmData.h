
  // DATA
  // very first word in the string - data block identifier 
  //        1 -   black stopper, before repairs 
  //        2 -   black stopper, after repairs 
  //        3 -   clear exit (no black stopper)
  // Clear side: lines 1-11
  // line 12 - Y-positions
  // Fiber side : lines 13-23
  // col 0 - R(Y=3) =  meas(A)/mes(B)
  // col 1 - meas(B)
  // col 2 - X (1-11, step 1)
  // col 3 - meas(Y=3)
  // col 4 - meas(Y=5)
  // col .....
  // col 11- meas(Y=19) - A
  // col 12- meas(Y=19) - B
  // col 13- R(Y=19) =  meas(A)/mes(B)
  // col 14-16 Comments

//  tile with black stopper before repairs
TString rd1 = "1,B,,A,,,,,,,,A,B,,,,,,,
,,11,,,,,,502,506,520,360,341,0.947222222,,,,,,
1.057057057,352,10,333,395,435,468,468,550,600,566,458,450,0.982532751,,LED,,,,
1.021108179,387,9,379,464,618,600,725,597,604,716,737,816,1.107191316,,Amplitude,4.1,,,
0.940758294,397,8,422,608,579,572,545,506,485,633,483,479,0.991718427,,width,140ns,,,
1.1003861,570,7,518,566,487,464,441,475,470,712,533,562,1.054409006,,,,,,
1.075,645,6,600,445,447,429,412,443,422,775,1810,1860,1.027624309,,,,64141,92,697.1847826
1.025270758,568,5,554,456,412,452,460,402,479,1830,1980,1980,1,,(Saturated),,,,
1.063559322,502,4,472,579,472,485,445,525,745,1980,1760,1720,0.977272727,,,,,,
1.092896175,400,3,366,581,608,695,691,762,1007,1510,1180,1050,0.889830508,,,,,,
1.055045872,345,2,327,437,464,522,489,510,641,954,1250,1200,0.96,,,,,,
1.032258065,320,1,310,404,445,493,418,422,447,595,870,854,0.981609195,,,,,,
,,,3,5,7,9,11,13,15,17,19,,,,,,,,
0.97382199,186,1,191,222,227,258,213,197,175,130,89,91,1.02247191,,,,,,
0.985645933,206,2,209,263,294,269,238,232,206,153,91,90.6,0.995604396,,LED,,,,
0.959641256,214,3,223,1050,1970,1970,1980,1970,1960,177,103,100,0.970873786,,Amplitude,3.9,,,
0.977653631,350,4,358,325,308,284,280,255,244,266,112,111,0.991071429,,width,150ns,,,
0.885321101,579,5,654,289,256,276,258,246,231,1980,140,142,1.014285714,,,,,,
0.695412844,758,6,1090,288,289,271,282,251,251,1930,166,166,1,,,,,,
0.786086957,452,7,575,327,284,282,262,263,252,1510,215,212,0.986046512,,(Saturated),,,,
0.854237288,252,8,295,1000,345,336,304,244,333,377,242,236,0.975206612,,,,,,
0.924170616,195,9,211,525,292,360,1990,1960,373,445,1320,1420,1.075757576,,,,,,
0.892156863,182,10,204,221,252,266,283,298,257,360,243,258,1.061728395,,,,,,
,,11,,,,,,233,258,256,180,177,0.983333333,,,,,,";

//  painted (repaired) tile with black stopper
TString rd2= "2,,,A,,,,,,,,A,,,,,
,,11,,,,,,,,,,,,,,
,,10,105,,,,,,,,32,,,,LED,
,,9,110,,,183,225,189,218,206,35,,,,Amplitude,3.45
,,8,128,132,182,196,181,166,155,188,125,,,,width,70ns
,,7,152,183,165,,141,,137,203,105,,,,,
,,6,175,205,,,,,,19.5,21.8,,,,,
,,5,151,176,,,,,,18.5,19.7,,,,(Saturated),
,,4,144,156,173,151,166,,21,21,19,,,,,
,,3,110,,178,194,191,242,24,19.2,21,,,,,
,,2,97,,,,160,44,21,17.5,18.8,,,,,
,,1,89,,,,,,,18,18.6,,,,,
,,,3,5,7,9,11,13,15,17,19,,,,,
,,1,76,,,,,,79,18,60,,,,,
,,2,87,,,,,,94,76,63,,,,,
,,3,92,22,40,25,19.7,21,25,82,68,,,,,
,,4,131,125,,,,,115,17,83,,,,,
,,5,64,125,,,,,107,22,92,,,,,
,,6,21.1,37,,,,,116,20.5,91,,,,,
,,7,162,22,,,,,,31.5,111,,,,,
,,8,115,21,68,137,,,,135,112,,,,,
,,9,91,200,112,28,23.4,63,21,35,27,,,,,
,,10,80,,,,,110,106,120,108,,,,,
,,11,,,,,,,,,93,,,,,";


// tile without black stopper
TString rd3 = "3,B,,A,,,,,,,,A,B,,,,,
,,11,,,,,147,165,180,232,242,208,1.163461538,,LED,,
1.0625,96,10,102,130,153,173,191,201,237,269,337,300,1.123333333,,,A,B
1.058333333,120,9,127,163,211,233,262,322,335,391,439,385,1.14025974,,Amplitude[V],4.1,1.83
1.072992701,137,8,147,226,231,198,198,209,239,276,294,253,1.162055336,,width[ns],140,150
1.058823529,170,7,180,181,167,165,161,157,208,252,213,188,1.132978723,,,,
1.197674419,172,6,206,172,154,152,154,145,159,245,179,153,1.169934641,,,,
1.056179775,178,5,188,162,148,137,150,156,155,216,158,130,1.215384615,,(Saturated),,
1.096774194,155,4,170,168,151,162,153,172,183,201,141,125,1.128,,,,
1.103448276,116,3,128,209,220,215,188,198,194,161,118,100,1.18,,,,
1.052173913,115,2,121,142,179,180,175,176,156,129,119,100,1.19,,,,
1.136842105,95,1,108,135,142,145,140,138,122,118,98,84,1.166666667,,,,
,,,3,5,7,9,11,13,15,17,19,,,,,,
1.112244898,49,1,54.5,101,75,60,68,65,59.5,47.5,35,36,0.972222222,,LED,,
1.19047619,52.5,2,62.5,94,87,88,84,83,85,61,44,47,0.936170213,,,A,B
1.19205298,60.4,3,72,1970,166,177,195,258,983,74,48.3,46.2,1.045454545,,Amplitude[V],3.9,1.8
1.114882507,76.6,4,85.4,100,84,83,88,95,93,131,61.2,53.3,1.148217636,,width[ns],150,150
1.079268293,82,5,88.5,87,80,70,66,68,76.6,720,79,80,0.9875,,,,
1.111111111,117,6,130,88,80,67,74,76,86.4,161,89.5,78,1.147435897,,,,
1.06779661,88.5,7,94.5,98.7,98,80,88,83,87.5,209,126,111,1.135135135,,(Saturated),,
1.019553073,71.6,8,73,333,132,111,112,111,113,152,145,135,1.074074074,,,,
1.211498973,48.7,9,59,94,110,112,168,800,1930,221,1930,712,2.710674157,,,,
1.157894737,47.5,10,55,64,74.5,73,97,107,123,140,184,182,1.010989011,,,,
,,11,,,,,75,76,98,105,115,112,1.026785714,,,,";


//  first measurement made on a tile with the stopper (or formatting problems)
TString rd = ",B,,A,,,,,,,,A,B,,,,
,,11,,,,,,502,506,520,360,341,0.947222222,,,
1.057057057,352,10,333,395,435,468,468,550,600,566,458,450,0.982532751,,LED,
1.021108179,387,9,379,464,618,600,725,597,604,716,737,816,1.107191316,,Amplitude,4.1,
0.940758294,397,8,422,608,579,572,545,506,485,633,483,479,0.991718427,,width,ns140,
1.1003861,570,7,518,566,487,464,441,475,470,712,533,562,1.054409006,,,
1.075,645,6,600,445,447,429,412,443,422,775,1810,1860,1.027624309,,,
1.025270758,568,5,554,456,412,452,460,402,479,1830,1980,1980,1,,Saturated,
1.063559322,502,4,472,579,472,485,445,525,745,1980,1760,1720,0.977272727,,,
1.092896175,400,3,366,581,608,695,691,762,1007,1510,1180,1050,0.889830508,,,
1.055045872,345,2,327,437,464,522,489,510,641,954,1250,1200,0.96,,,
1.032258065,320,1,310,404,445,493,418,422,447,595,870,854,0.981609195,,,
,,,3,5,7,9,11,13,15,17,19,,,,,
0.97382199,186,1,191,222,227,258,213,197,175,130,89,91,1.02247191,,,
0.985645933,206,2,209,263,294,269,238,232,206,153,91,90.6,0.995604396,,LED,
0.959641256,214,3,223,1050,1970,1970,1980,1970,1960,177,103,100,0.970873786,,Amplitude,3.9,
0.977653631,350,4,358,325,308,284,280,255,244,266,112,111,0.991071429,,width,ns150,
0.885321101,579,5,654,289,256,276,258,246,231,1980,140,142,1.014285714,,,
0.695412844,758,6,1090,288,289,271,282,251,251,1930,166,166,1,,,
0.786086957,452,7,575,327,284,282,262,263,252,1510,215,212,0.986046512,,Saturated,
0.854237288,252,8,295,1000,345,336,304,244,333,377,242,236,0.975206612,,,
0.924170616,195,9,211,525,292,360,1990,1960,373,445,1320,1420,1.075757576,,,
0.892156863,182,10,204,221,252,266,283,298,257,360,243,258,1.061728395,,,
  ,,11,,,,,,233,258,256,180,177,0.983333333,,,";

Double_t tScale[] = {1.44, 1., 1.};
Double_t exitX = 9;
Double_t exitY = 9;
// list of positions with visible effect of fiber
Int_t fib[] = {
  0,0,0,1,1,1,1,1,0,0,0,   0,0,1,0,0,0,0,1,0,0,0,   0,0,1,0,0,0,0,0,1,0,0,
  0,0,1,0,0,0,0,0,1,0,0,   0,0,1,0,0,0,0,0,1,0,0,   0,0,1,0,0,0,0,0,1,0,0,
  0,0,1,0,0,0,0,0,1,0,0,   0,0,1,1,1,1,1,1,1,0,0,   0,0,0,0,0,0,0,0,1,0,0};

// list of positions with coating subject to damage
Int_t dam[] = {
  0,0,0,0,0,0,0,0,0,0,1,   0,0,0,0,0,0,0,0,0,0,1,   0,0,0,0,0,0,0,0,0,0,1,
  0,0,0,0,0,0,0,0,0,0,1,   0,0,0,0,0,0,0,0,0,0,1,   0,0,0,0,0,0,0,0,0,0,0,
  0,1,1,1,0,0,0,0,0,0,0,   1,1,1,1,1,1,1,0,0,0,0,   1,1,1,1,1,1,1,0,0,0,0};

 