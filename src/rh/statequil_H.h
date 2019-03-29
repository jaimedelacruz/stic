#ifndef STATEH_H
#define STATEH_H

void CollisionRateOne(struct Atom *atom, char **fp_atom, int k);
void statEquil_H(Atom *atom, int isum, int mali_iter);
void SetLTEQuantitiesOne(Atom *atom, int k);
void getfjk2(Element *element, double ne, int k, double *fjk, double *dfjk);
double getKuruczpf2(Element *element, int stage, int k);
void FixedRateOne(Atom *atom, int k);


#endif
