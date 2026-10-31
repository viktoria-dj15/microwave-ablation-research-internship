// R35 + lokalno usitnjavanje koaksijalnog porta 201 (hPort=0.10 mm).
// Referenca: antena_V3_2_tumor107_R35_scaling_fix.geo.
// Ne menjati geometriju i ostale meshing parametre za uporedni test.
SetFactory("OpenCASCADE");
// KRITICNO: prilozeni BREP je VEC U METRIMA!
// Stari tumor_h1.5_m(2).geo postavlja OCCScaling=0.001.
// Ako Gmsh zapamti opciju, ShapeFromFile ponovo smanji tumor 1000x.
Geometry.OCCScaling = 1;

// DIMENZIJE 
mm=0.001;
rUnutrasnji = 0.25*mm;
rIzolator = 0.8*mm;
rSpoljasnji = 1*mm;
duzina = 50*mm;

G1 = 1*mm;
G2 = 4*mm;
L = 8*mm;

duzinaKapice = 0.5*mm;
duzinaVrha = 1*mm;

// KATETER
debljinaKatetera = 0.2*mm;
rKateter = rSpoljasnji + debljinaKatetera;

// OKOLNO TKIVO
// Probni R=35 mm: ~15 mm minimalnog radijalnog prostora od tumora do ABC.
// Studija osetljivosti: uporediti sa prethodnom mrezom R=27 mm.
rDomen = 35*mm;
produzetakDomena = 20*mm;

// Krajnji polozaji
xKraj = duzina + duzinaKapice + duzinaVrha;
xKrajKatetera = xKraj + debljinaKatetera;
xKrajDomena = xKrajKatetera + produzetakDomena;

eps = 1e-5;

// POLOZAJI PROREZA — racunaju se od kraja antene
xG2 = duzina - G2;
xG1 = xG2 - L - G1;


// UNUTRASNJI PROVODNIK
unutrasnji = newv;
Cylinder(unutrasnji) = {0, 0, 0, duzina, 0, 0, rUnutrasnji, 2*Pi};

// IZOLATOR
izolatorPun = newv;
Cylinder(izolatorPun) = {0, 0, 0, duzina, 0, 0, rIzolator, 2*Pi};

izolator[] = BooleanDifference{Volume{izolatorPun}; Delete;} {Volume{unutrasnji}; };


// SPOLJASNJI METALNI PROVODNIK
spoljasnjiPun = newv;
Cylinder(spoljasnjiPun) = {0, 0, 0, duzina, 0, 0, rSpoljasnji, 2*Pi};

pomocni = newv;
Cylinder(pomocni) = {0, 0, 0, duzina, 0, 0, rIzolator, 2*Pi};

cev[] = BooleanDifference{ Volume{spoljasnjiPun}; Delete;}{Volume{pomocni}; Delete;};



// DVA PROREZA
rez1 = newv;
Cylinder(rez1) = {xG1, 0, 0, G1, 0, 0, 1.1*rSpoljasnji, 2*Pi};

rez2 = newv;
Cylinder(rez2) = {xG2, 0, 0, G2, 0, 0, 1.1*rSpoljasnji, 2*Pi};

spoljasnji[] = BooleanDifference{ Volume{cev[]}; Delete; }{ Volume{rez1, rez2}; Delete;};


// ZAPREMINA PRVOG SLOTA
slot1Pun = newv;
Cylinder(slot1Pun) = {xG1, 0, 0, G1, 0, 0, rSpoljasnji, 2*Pi};

slot1Rupa = newv;
Cylinder(slot1Rupa) = {xG1, 0, 0, G1, 0, 0, rIzolator, 2*Pi};

slot1[] = BooleanDifference{ Volume{slot1Pun}; Delete;}{ Volume{slot1Rupa}; Delete;};


// ZAPREMINA DRUGOG SLOTA
slot2Pun = newv;
Cylinder(slot2Pun) = { xG2, 0, 0, G2, 0, 0, rSpoljasnji, 2*Pi};

slot2Rupa = newv;
Cylinder(slot2Rupa) = { xG2, 0, 0, G2, 0, 0, rIzolator, 2*Pi};

slot2[] = BooleanDifference{ Volume{slot2Pun}; Delete; }{ Volume{slot2Rupa}; Delete; };


// KAPICA I VRH
kapica = newv;
Cylinder(kapica) = {
  duzina, 0, 0,
  duzinaKapice, 0, 0,
  rSpoljasnji, 2*Pi
};

vrh = newv;
Cone(vrh) = {
  duzina + duzinaKapice, 0, 0,
  duzinaVrha, 0, 0,
  rSpoljasnji, 0, 2*Pi
};

provodnikSaKapicom[] = BooleanUnion{ Volume{unutrasnji}; Delete; }{ Volume{kapica, vrh}; Delete; };

// KATETER: SPOLJASNJI OBLIK

kateterCilindar = newv;
Cylinder(kateterCilindar) = {
  0, 0, 0,
  duzina + duzinaKapice, 0, 0,
  rKateter, 2*Pi
};

kateterVrh = newv;
Cone(kateterVrh) = {
  duzina + duzinaKapice, 0, 0,
  duzinaVrha + debljinaKatetera, 0, 0,
  rKateter, 0, 2*Pi
};

kateterPun[] = BooleanUnion{
  Volume{kateterCilindar}; Delete;
}{
  Volume{kateterVrh}; Delete;
};

// Iz spoljasnjeg oblika uklanjamo prostor koji zauzima antena.
// Antenu zadrzavamo: u drugoj zagradi nema Delete.
kateter[] = BooleanDifference{
  Volume{kateterPun[]}; Delete;
}{
  Volume{
    provodnikSaKapicom[],
    spoljasnji[],
    izolator[],
    slot1[],
    slot2[]
  };
};

// OKOLNI RACUNARSKI DOMEN

domenPun = newv;
Cylinder(domenPun) = {
  0, 0, 0,
  xKrajDomena, 0, 0,
  rDomen, 2*Pi
};

// Tkivo zauzima prostor izvan antene i katetera.
// Antenu i kateter zadrzavamo.
tkivo[] = BooleanDifference{
  Volume{domenPun}; Delete;
}{
  Volume{
    provodnikSaKapicom[],
    spoljasnji[],
    izolator[],
    slot1[],
    slot2[],
    kateter[]
  };
};

// Prvo uskladjujemo originalnu antenu, kateter i tkivo BEZ tumora.
// Originalna V3.2 geometrija ostaje nepromenjena, osim R=35 mm.
Coherence;

// Originalni tumor h=1.5 (3580 ravnih CAD lica) prethodno je
// transliran na x=43.5 mm i iz njega je ODUZET SPOLJASNJI OMOTAC
// katetera (cilindar r=1.2 mm + konusni vrh).
// Tako tumor ne moze preklapati antenu/kateter.
// Ovaj BREP JE POTREBAN u istom folderu kao i ovaj .geo.
// Reset neposredno pre uvoza (zastita od prethodno ucitanih .geo fajlova)
Geometry.OCCScaling = 1;
vTumorUlaz[] = ShapeFromFile("tumor_h1_5_107_izvan_katetera_m.brep");
If (#vTumorUlaz[] != 1)
  Error("Ocekuje se jedan zatvoren CAD solid tumora. Proveri BREP fajl!");
EndIf
bbTumor() = BoundingBox Volume{vTumorUlaz[]};
Printf("TUMOR 107 BREP: X[%g,%g] mm, Y[%g,%g] mm, Z[%g,%g] mm", 
  bbTumor(0)/mm, bbTumor(3)/mm,
  bbTumor(1)/mm, bbTumor(4)/mm,
  bbTumor(2)/mm, bbTumor(5)/mm);
// Ocekivani X: priblizno 27.7..59.3 mm (nas tumor h=1.5).
// STOP pre BooleanFragments ako je doslo do ponovnog skaliranja!
If (bbTumor(0) < 20*mm)
  Error("POGRESNO SKALIRANJE TUMORA: Xmin manji od 20 mm. Zatvori Gmsh, ponovo otvori ovaj .geo; OCCScaling mora biti 1!");
EndIf
If (bbTumor(3) > 65*mm)
  Error("NEOCEKIVAN POLOZAJ TUMORA: Xmax veci od 65 mm. Prekinuti pre meshinga!");
EndIf
If ((bbTumor(3)-bbTumor(0)) < 25*mm)
  Error("POGRESNA VELICINA TUMORA: sirina po X manja od 25 mm!");
EndIf
Printf("PROVERA RAZMERA TUMORA: OK, OCCScaling=1, X dimenzije u mm proverene");

// SVE ORIGINALNE ZAPREMINE, osim uvezenog tumora, dele se
// JEDNOM OCC OPERACIJOM sa tumorom: zajednicki interfejsi 106/107
// i 105/107 postaju konformni, nema preklopljenih tetraedara.
vPostojece[] = Volume{:};
vPostojece[] -= vTumorUlaz[];
vSveFragmenti[] = BooleanFragments{
  Volume{vPostojece[]}; Delete;
}{
  Volume{vTumorUlaz[]}; Delete;
};
Printf("KONFORMNA PODELA: %g CAD fragment zapremina", #vSveFragmenti[]);
// VAŽNO: posle ovog koraka NE pozivati Coherence ponovo:
// svi odgovarajuci tagovi se traze ponovo ispod.


// SVE ZAPREMINE ANTENE
xKraj = duzina + duzinaKapice + duzinaVrha;

vAntena[] = Volume In BoundingBox{
  -eps, -rSpoljasnji-eps, -rSpoljasnji-eps,
  xKraj+eps, rSpoljasnji+eps, rSpoljasnji+eps
};

// PRONALAZENJE PRVOG SLOTA
vSlot1[] = Volume In BoundingBox{
  xG1-eps, -rSpoljasnji-eps, -rSpoljasnji-eps,
  xG1+G1+eps, rSpoljasnji+eps, rSpoljasnji+eps
};

Physical Volume("Slot1", 103) = {vSlot1[]};

// PRONALAZENJE DRUGOG SLOTA
vSlot2[] = Volume In BoundingBox{
  xG2-eps, -rSpoljasnji-eps, -rSpoljasnji-eps,
  xG2+G2+eps, rSpoljasnji+eps, rSpoljasnji+eps
};

Physical Volume("Slot2", 104) = {vSlot2[]};

// IZOLATOR
vIzolator[] = Volume In BoundingBox{
  -eps, -rIzolator-eps, -rIzolator-eps,
  duzina+eps, rIzolator+eps, rIzolator+eps
};

Physical Volume("Izolator", 102) = {vIzolator[]};

// IZDVAJANJE METALA
vMetal[] = vAntena[];
vMetal[] -= vIzolator[];
vMetal[] -= vSlot1[];
vMetal[] -= vSlot2[];

Physical Volume("Metal", 101) = {vMetal[]};

// KATETER: kutija obuhvata antenu i kateter
vKateter[] = Volume In BoundingBox{
  -eps, -rKateter-eps, -rKateter-eps,
  xKrajKatetera+eps, rKateter+eps, rKateter+eps
};

// Uklanjamo antenu: ostaje kateter
vKateter[] -= vAntena[];

Physical Volume("Kateter", 105) = {vKateter[]};


// TUMOR: prethodni BREP bbox ostaje zapamcen PRE fragmentacije.
// Region 107 je geometrijski sadrzan u toj kutiji; izbacujemo
// eventualne unutrasnje delove antene/katetera.
// Uvecaj bbox za 0.15 mm zbog OCC tolerancija nakon fragmentacije.
// Ovo ne pomera tumor i ne menja njegovu CAD geometriju.
epsTumor = 0.15*mm;
vTumor107[] = Volume In BoundingBox{
  bbTumor(0)-epsTumor, bbTumor(1)-epsTumor, bbTumor(2)-epsTumor,
  bbTumor(3)+epsTumor, bbTumor(4)+epsTumor, bbTumor(5)+epsTumor
};
vTumor107[] -= vAntena[];
vTumor107[] -= vKateter[];
If (#vTumor107[] < 1)
  Error("TUMOR 107 NEDOSTAJE POSLE FRAGMENTACIJE - NE GENERISI MREZU. Posalji KONFORMNA PODELA i bbox ispis!");
EndIf
Physical Volume("Tumor", 107) = {vTumor107[]};
Printf("DIJAGNOSTIKA: tumor 107 ima %g CAD zapremina", #vTumor107[]);

// TKIVO: kutija obuhvata ceo model
vTkivo[] = Volume In BoundingBox{
  -eps, -rDomen-eps, -rDomen-eps,
  xKrajDomena+eps, rDomen+eps, rDomen+eps
};

// Uklanjamo antenu i kateter: ostaje tkivo
vTkivo[] -= vAntena[];
vTkivo[] -= vKateter[];
vTkivo[] -= vTumor107[];
If (#vTkivo[] < 1)
  Error("ZDRAVO TKIVO 106 NEDOSTAJE posle OCC fragmentacije!");
EndIf
Physical Volume("Tkivo", 106) = {vTkivo[]};
Printf("DIJAGNOSTIKA: zdravo tkivo 106 ima %g CAD zapremina", #vTkivo[]);



// PORT: ULAZNI PRSTEN IZOLATORA NA x = 0
// Povrsine na ulazu, do spoljasnjeg poluprecnika izolatora
sUlaz[] = Surface In BoundingBox{
  -eps, -rIzolator-eps, -rIzolator-eps,
   eps,  rIzolator+eps,  rIzolator+eps
};

// Centralni metalni disk na ulazu
sCentralniDisk[] = Surface In BoundingBox{
  -eps, -rUnutrasnji-eps, -rUnutrasnji-eps,
   eps,  rUnutrasnji+eps,  rUnutrasnji+eps
};


// Iz liste uklanjamo disk: ostaje prsten izolatora
sPort[] = sUlaz[];
sPort[] -= sCentralniDisk[];

Physical Surface("Port", 201) = {sPort[]};


// PEC: DODIR METALA SA IZOLATOROM I SLOTOVIMA

// Sve granicne povrsine metala
sGranicaMetala[] = Abs(Boundary{
  Volume{vMetal[]};
});

// Granicne povrsine izolatora, slotova i katetera
sGranicaDielektrika[] = Abs(Boundary{
  Volume{vIzolator[], vSlot1[], vSlot2[], vKateter[]};
});

// Povrsine dielektrika koje ne dodiruju metal
sNisuMetal[] = sGranicaDielektrika[];
sNisuMetal[] -= sGranicaMetala[];

// Zajednicke povrsine dielektrika i metala
sPEC[] = sGranicaDielektrika[];
sPEC[] -= sNisuMetal[];

Physical Surface("PEC", 202) = {sPEC[]};
// ABSORBING: SPOLJASNJA GRANICA TKIVA

// Sve granicne povrsine tkiva
sGranicaTkiva[] = Abs(Boundary{Volume{vTkivo[]};});

// Tumor/jetra su UNUTRASNJE granice: NIKAKO ne smeju biti ABC!
sGranicaTumora[] = Abs(Boundary{Volume{vTumor107[]};});

// Sve granicne povrsine katetera
sGranicaKatetera[] = Abs(Boundary{Volume{vKateter[]};});

// Uklanjamo dodir tkiva sa kateterom.
// Ostaju spoljasnje povrsine cilindra tkiva.
sAbsorbing[] = sGranicaTkiva[];
sAbsorbing[] -= sGranicaKatetera[];
sAbsorbing[] -= sGranicaTumora[];

// Iskljucujemo povrsine sa drugim granicnim uslovima
sAbsorbing[] -= sPort[];
sAbsorbing[] -= sPEC[];

Physical Surface("Absorbing", 203) = {sAbsorbing[]};



// V3.2: OPTIMIZOVANA MREZA ZA MPI + MUMPS
// Fizicka geometrija i Physical Volume/Surface grupe ostaju iste kao u V2.
// Sve dimenzije u ovom bloku su u metrima (mm = 0.001).
//
// Ideja:
//  - zadrzati 0.12 mm oko uskog unutrasnjeg provodnika (bolji kruzni presek),
//    ali znatno suziti oblast u kojoj se primenjuje ta vrednost;
//  - zadrzati 0.15 mm neposredno oko oba slota;
//  - opustiti mrezni korak u ostatku antene/katetera i prelaznoj zoni tkiva;
//  - ne usitnjavati udaljeno tkivo zbog prevelikog Thickness polja.
//
// V3 je probna optimizacija: pre MUMPS-a proveriti broj elemenata,
// fizicke grupe, kvalitet i poprecni presek centralnog provodnika.

hBlizu = 0.55*mm;       // V2: 0.40 mm (srednji/udaljeniji delovi antene)
hDaleko = 4.0*mm;      // ne menjamo spoljni domen
hUnutrasnji = 0.12*mm; // V2: 0.12 mm - cuvamo finu konturu provodnika
hSlot = 0.15*mm;       // V2: 0.15 mm - cuvamo rezoluciju oba slota
hPort = 0.10*mm;       // NOVO: lokalna velicina elemenata oko porta 201 na x=0

// Field 1: osnovni sloj uz kateter, sa kracim prelazom u okolno tkivo.
// V2 je imao Thickness = 5 mm i dodatno usitnjavao znatnu oblast tkiva.
marginaMreze = 0.15*mm;
Field[1] = Box;
Field[1].VIn = hBlizu;
Field[1].VOut = hDaleko;
Field[1].XMin = -marginaMreze;
Field[1].XMax = xKrajKatetera + marginaMreze;
Field[1].YMin = -rKateter - marginaMreze;
Field[1].YMax =  rKateter + marginaMreze;
Field[1].ZMin = -rKateter - marginaMreze;
Field[1].ZMax =  rKateter + marginaMreze;
Field[1].Thickness = 1.2*mm;

// Field 2: veoma uska kutija oko centralnog metalnog provodnika.
// U V2 je obuhvat bio +/-0.40 mm po y i z, uz Thickness 0.30 mm.
// Sada je obuhvat +/-0.30 mm, uz kraci prelaz od 0.20 mm.
// Time 0.12 mm ostaje oko granice provodnika (r=0.25 mm),
// a ne namece se nepotrebno kroz veliki deo izolatora.
rZonaUnutrasnji = rUnutrasnji + 0.05*mm;
Field[2] = Box;
Field[2].VIn = hUnutrasnji;
Field[2].VOut = hDaleko;
Field[2].XMin = -0.10*mm;
Field[2].XMax = duzina + 0.10*mm;
Field[2].YMin = -rZonaUnutrasnji;
Field[2].YMax =  rZonaUnutrasnji;
Field[2].ZMin = -rZonaUnutrasnji;
Field[2].ZMax =  rZonaUnutrasnji;
Field[2].Thickness = 0.20*mm;

// Field 3: prvi slot, kratak prorez x = 37..38 mm.
// V2: aksijalna margina 0.40 mm i radijalna kutija +/-1.40 mm.
// V3 suzava samo zonu maksimalnog usitnjavanja, ne sam slot.
rZonaSlot = rKateter + 0.05*mm;
Field[3] = Box;
Field[3].VIn = hSlot;
Field[3].VOut = hDaleko;
Field[3].XMin = xG1 - 0.25*mm;
Field[3].XMax = xG1 + G1 + 0.25*mm;
Field[3].YMin = -rZonaSlot;
Field[3].YMax =  rZonaSlot;
Field[3].ZMin = -rZonaSlot;
Field[3].ZMax =  rZonaSlot;
Field[3].Thickness = 0.45*mm;

// Field 4: drugi slot x = 46..50 mm, kapica i metalni vrh.
Field[4] = Box;
Field[4].VIn = hSlot;
Field[4].VOut = hDaleko;
Field[4].XMin = xG2 - 0.25*mm;
Field[4].XMax = xKraj + 0.15*mm;
Field[4].YMin = -rZonaSlot;
Field[4].YMax =  rZonaSlot;
Field[4].ZMin = -rZonaSlot;
Field[4].ZMax =  rZonaSlot;
Field[4].Thickness = 0.45*mm;

// V3.1: Dodatno usitnjavanje samo u tankom prstenu uz SPOLJASNJU
// granicu izolatora (r = 0.80 mm), duz x = 0..50 mm.
// Ovo ne menja geometriju, vec samo ciljnu velicinu elemenata.
// Na samoj granici r=0.80 mm cilj je 0.24 mm.
// Van uskog radijalnog pojasa polje vraca hDaleko, tj. ne utice na mrezu.
// Funkcije Sqrt, Min i Step su podrzane u Gmsh MathEval izrazu.
hIzolatorSpolja = 0.24*mm;
sirinaPrstena = 0.18*mm;

Field[6] = MathEval;
Field[6].F = Sprintf(
  "%g-(%g-%g)*(1-Min(1,((Sqrt(y*y+z*z)-%g)/%g)^2))*Step(x)*Step(%g-x)",
  hDaleko, hDaleko, hIzolatorSpolja, rIzolator, sirinaPrstena, duzina
);

// V3.2: dodatno fino usitnjavanje na dve SPOLJASNJE granice.
// r = 1.00 mm: spoljasnji metal/kateter (vazan i za PEC povrsine).
// r = 1.20 mm: kateter/tkivo (ukupni kruzni presek ranije -2.69%).
// Kroz metale i izolator NE menjamo CAD geometriju ni fizicke oznake.
hSpoljasnjiMetal = 0.25*mm;
hSpoljasnjiKateter = 0.30*mm;
sirinaSpoljasnjihPrstenova = 0.18*mm;

// Fino polje uz r=1.00 mm, na cilindricnom delu antene 0..50 mm.
Field[7] = MathEval;
Field[7].F = Sprintf(
  "%g-(%g-%g)*(1-Min(1,((Sqrt(y*y+z*z)-%g)/%g)^2))*Step(x+%g)*Step(%g-x)",
  hDaleko, hDaleko, hSpoljasnjiMetal,
  rSpoljasnji, sirinaSpoljasnjihPrstenova,
  0.10*mm, duzina+0.10*mm
);

// Fino polje uz r=1.20 mm, na cilindricnom delu katetera 0..50.5 mm.
Field[8] = MathEval;
Field[8].F = Sprintf(
  "%g-(%g-%g)*(1-Min(1,((Sqrt(y*y+z*z)-%g)/%g)^2))*Step(x+%g)*Step(%g-x)",
  hDaleko, hDaleko, hSpoljasnjiKateter,
  rKateter, sirinaSpoljasnjihPrstenova,
  0.10*mm, duzina+duzinaKapice+0.10*mm
);

// Min uzima najfiniju ciljnu velicinu svih lokalnih polja.
// TUMOR h=1.5: originalna CAD povrsina ostaje nepromenjena;
// gustina novih tetraedara kontrolise se oko njenih faceta.
// Velicina 1.4 mm blizu tumora, postepeno 4 mm daleko od njega.
Field[9] = Distance;
Field[9].SurfacesList = {sGranicaTumora[]};
Field[9].Sampling = 3;
Field[10] = Threshold;
Field[10].InField = 9;
Field[10].SizeMin = 1.4*mm;
Field[10].SizeMax = hDaleko;
Field[10].DistMin = 0.3*mm;
Field[10].DistMax = 3.5*mm;

// NOVO: samo lokalno usitnjavanje ulaznog prstena koaksijalnog porta 201.
// Port je u ravni x=0, r=0.25..0.80 mm. Box obuhvata ceo
// prsten i samo prvih 0.8 mm duz koaksijalnog voda.
// Cilj je povecati broj trouglova 201 (R35 referenca: 120),
// bez globalnog usitnjavanja tumora ili spoljasnjeg tkiva.
// Ovo polje ne menja CAD, port, grupe ni granične uslove.
Field[11] = Box;
Field[11].VIn = hPort;
Field[11].VOut = hDaleko;
Field[11].XMin = -0.10*mm;
Field[11].XMax = 0.80*mm;
Field[11].YMin = -rIzolator - 0.05*mm;
Field[11].YMax =  rIzolator + 0.05*mm;
Field[11].ZMin = -rIzolator - 0.05*mm;
Field[11].ZMax =  rIzolator + 0.05*mm;
Field[11].Thickness = 0.25*mm;

Field[5] = Min;
Field[5].FieldsList = {1, 2, 3, 4, 6, 7, 8, 10, 11};
Background Field = 5;

// Iskljucujemo ostale automatske izvore velicine mreze,
// da nam postavljeni Box/Min profil ostane kontrolisan.
Mesh.MeshSizeFromPoints = 0;
Mesh.MeshSizeFromCurvature = 0;
Mesh.MeshSizeExtendFromBoundary = 0;

// Ovaj parametar zadrzavamo kao u V2 radi aproksimacije malih krugova.
Mesh.MinimumCirclePoints = 24;

Mesh.MeshSizeMin = hPort; // NOVO: omogucava korak 0.10 mm kod porta 201
Mesh.MeshSizeMax = hDaleko;

Printf("Slot1: %g zapremina", #vSlot1[]);
Printf("Slot2: %g zapremina", #vSlot2[]);
Printf("Izolator: %g zapremina", #vIzolator[]);
Printf("Antena ukupno: %g zapremina", #vAntena[]);
Printf("Metal: %g zapremina", #vMetal[]);
Printf("Ulaz: %g povrsina", #sUlaz[]);
Printf("Centralni disk: %g povrsina", #sCentralniDisk[]);
Printf("Port: %g povrsina", #sPort[]);
Printf("PORT REFINEMENT: hPort=%g mm, x-zona [-0.10, 0.80] mm; R=35 mm", hPort/mm);
Printf("PEC: %g povrsina", #sPEC[]);
Printf("Kateter: %g zapremina", #vKateter[]);
Printf("Tkivo: %g zapremina", #vTkivo[]);
Printf("Absorbing: %g povrsina", #sAbsorbing[]);

vNepokriveno[] = Volume{:};
vNepokriveno[] -= vMetal[];
vNepokriveno[] -= vIzolator[];
vNepokriveno[] -= vSlot1[];
vNepokriveno[] -= vSlot2[];
vNepokriveno[] -= vKateter[];
vNepokriveno[] -= vTkivo[];
vNepokriveno[] -= vTumor107[];

Printf("Zapremine bez Physical Volume grupe: %g", #vNepokriveno[]);

Show "*";


// Dodatne dijagnostike (sam CAD, PRE meshing):
Printf("DIJAGNOSTIKA PRE MESH: Tumor 107: %g volume, jetra 106: %g volume", #vTumor107[], #vTkivo[]);
Printf("DIJAGNOSTIKA PRE MESH: Port201=%g PEC202=%g ABC203=%g povrsina", #sPort[], #sPEC[], #sAbsorbing[]);
vNeoznaceno[] = Volume{:};
vNeoznaceno[] -= vMetal[];
vNeoznaceno[] -= vIzolator[];
vNeoznaceno[] -= vSlot1[];
vNeoznaceno[] -= vSlot2[];
vNeoznaceno[] -= vKateter[];
vNeoznaceno[] -= vTkivo[];
vNeoznaceno[] -= vTumor107[];
If (#vNeoznaceno[] != 0)
  Printf("GRESKA: %g neoznacenih CAD zapremina", #vNeoznaceno[]);
  Error("Neki CAD regioni nemaju fizicku grupu");
EndIf
// Ne pustati MUMPS ako Gmsh javi 'No elements in volume'!
// Podesavanje Mesh.Algorithm3D po potrebi iz Gmsh GUI (HXT kao drugi pokusaj).
