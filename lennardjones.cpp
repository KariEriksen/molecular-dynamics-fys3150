#include "lennardjones.h"
#include "system.h"
#include <cmath>

double LennardJones::potentialEnergy() const
{

    return m_potentialEnergy;
}

double LennardJones::sigma() const
{
    return m_sigma;
}

void LennardJones::setSigma(double sigma)
{
    m_sigma = sigma;
}

double LennardJones::epsilon() const
{
    return m_epsilon;
}

void LennardJones::setEpsilon(double epsilon)
{
    m_epsilon = epsilon;
}

void LennardJones::calculateForces(System &system)
{
    m_potentialEnergy = 0; // Remember to compute this in the loop
    //Xforce = atom->position(xi - xj)*24*epsilon/rij^(2)*((sigma/rij)^6 - 2*(sigma/rij)^12);
    //System system;

    for (int i=0; i < system.atoms().size(); i++){
        for (int j=i+1; j < system.atoms().size(); j++) {

            Atom* atomi = system.atoms().at(i);
            Atom* atomj = system.atoms().at(j);
            vec3 dr = atomi->position - atomj->position;

            double dx = dr.x();
            double dy = dr.y();
            double dz = dr.z();

            double xSize = system.systemSize().x();
            double ySize = system.systemSize().y();
            double zSize = system.systemSize().z();

            if (dx > xSize/2) dx = dx - xSize;
            if (dy > ySize/2) dy = dy - ySize;
            if (dz > zSize/2) dz = dz - zSize;

            if (dx <= -xSize/2) dx = dx + xSize;
            if (dy <= -ySize/2) dy = dy + ySize;
            if (dz <= -zSize/2) dz = dz + zSize;

            double s = sigma();
            double sigma6 = s*s*s*s*s*s;
            double sigma12 = sigma6*sigma6;
            double rij2 = (dx*dx + dy*dy + dz*dz);
            double rij = sqrt(rij2);
            double rij6 = rij2*rij2*rij2;
            double rij12 = rij6*rij6;

            m_potentialEnergy += 4*epsilon()*((sigma12/rij12) - (sigma6/rij6));

            double Xforce = (dx/rij2)*24*epsilon()*((sigma6/rij6) - 2*(sigma12/rij12));
            double Yforce = (dy/rij2)*24*epsilon()*((sigma6/rij6) - 2*(sigma12/rij12));
            double Zforce = (dz/rij2)*24*epsilon()*((sigma6/rij6) - 2*(sigma12/rij12));

            vec3 dforce(Xforce, Yforce, Zforce);
            atomi->force -= dforce;
            atomj->force += dforce;



        }
    }
}
