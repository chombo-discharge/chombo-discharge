#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#include "Stencil.H"

#include "NamespaceHeader.H"

// IndexML:operator<<
ostream&
operator<< (ostream&       os,
            const IndexML& p)
{
  os << '(' << p.m_iv << ", lev:" << p.level() << /* ", b:" << p.block() << */ ')';
  if (os.fail())
    MayDay::Error("operator<<(ostream&,IndexML&) failed");
  return os;
}
bool
IndexML::operator> (const IndexML& p) const
{
  if ( m_lev < p.m_lev ) return true;
  else if ( m_lev > p.m_lev ) return false;
 
  if ( m_blockID < p.m_blockID ) return true;
  else if ( m_blockID > p.m_blockID ) return false;

  if ( m_iv[0] > p.m_iv[0] ) return true;
  else if (m_iv[0] < p.m_iv[0] ) return false;
#if CH_SPACEDIM > 1
  if ( m_iv[1] > p.m_iv[1] ) return true;
  else if (m_iv[1] < p.m_iv[1] ) return false;
#if CH_SPACEDIM > 2
  if ( m_iv[2] > p.m_iv[2] ) return true;
  else return false;
#endif
#endif
  return false;
}

bool
IndexML::operator< (const IndexML& p) const
{
  if ( m_lev > p.m_lev ) return true;
  else if ( m_lev < p.m_lev ) return false;

  if ( m_blockID > p.m_blockID ) return true;
  else if ( m_blockID < p.m_blockID ) return false;

  if ( m_iv[0] < p.m_iv[0] ) return true;
  else if (m_iv[0] > p.m_iv[0] ) return false;
#if CH_SPACEDIM > 1
  if ( m_iv[1] < p.m_iv[1] ) return true;
  else if (m_iv[1] > p.m_iv[1] ) return false;
#if CH_SPACEDIM > 2
  if ( m_iv[2] < p.m_iv[2] ) return true;
  else return false;
#endif
#endif
  
  return false;
}

bool
IndexML::operator==(const IndexML& p) const
{
  return (p.m_iv==m_iv && p.m_lev==m_lev && p.m_blockID==m_blockID);
}

bool
IndexML::operator!= (const IndexML& p) const
{
  return (p.m_iv!=m_iv || p.m_lev!=m_lev || p.m_blockID!=m_blockID);
}

// StencilTensorValue:operator<<
ostream&
operator<< (ostream&       os,
            const StencilTensorValue& p)
{
  if (p.m_dof==1)
    {
      os << p.value(0,0);
    }
  else
    {
      for (int i=0;i<p.m_dof;i++) 
        {
          for (int j=0;j<p.m_dof-1;j++) os << p.value(i,j) << ',';
          os << p.value(i,p.m_dof-1) << std::endl;
        }
    }
  if (os.fail())
    MayDay::Error("operator<<(ostream&,StencilTensorValue&) failed");
  return os;
}

// // StencilScalarValue::operator=
// StencilScalarValue& StencilScalarValue::operator=(const StencilScalarValue& p)
// {
//   m_val = p.m_val;
//   return *this;
// }

// // StencilScalarValue:operator<<
// ostream&
// operator<< (ostream&       os,
//             const StencilScalarValue& p)
// {
//   os << p.value();
//   if (os.fail())
//     MayDay::Error("operator<<(ostream&,StencilScalarValue&) failed");
//   return os;
// }

// A specialized operator that "distributes" a_sten[jiv] with a_new: a_sten += a_new*a_sten[jiv]. 
// Would like to remove a_sten[jiv] but that messes up iterators.
void
StencilProject(IndexML a_mliv, Vector<StencilNode> &a_scales, StencilTensor &a_sten)
{
  //CH_TIME("PetscCompGrid::projectStencil");
  StencilTensorValue ghostNode = a_sten[a_mliv]; // node getting deleted (distributed)
  // would like to remove the root but this messes up interators
  StencilTensor::iterator root = a_sten.find(a_mliv);
  // add scaled
  for (int i=0;i<a_scales.size();i++) 
    {
      StencilNode &target = a_scales[i];
      CH_assert(target.first != a_mliv); // this OK in theory but never done now
      StencilTensorValue &val = a_sten[target.first]; val.define(root->second); // only place where nodes are added
      val.apply(target.second,root->second);
    }
}

// axpy for tensor nodes: this = this + scale * node
void 
StencilTensorValue::apply(StencilTensorValue &a_scale, StencilTensorValue &a_node) 
{
  CHECK_DOF; // little DGEMM
  if (a_node.m_dof != m_dof) MayDay::Abort("StencilTensorValue::apply node DOF mismatch");
  if (a_node.m_dof != a_scale.m_dof && a_scale.m_dof != 1) MayDay::Abort("StencilTensorValue::apply scale DOF mismatch");
  if (a_scale.m_dof == 1) {
    for (int i=0;i<m_dof*m_dof;i++) m_val[i] += a_scale.m_val[0]*a_node.m_val[i];
  } else {
    for (int i=0;i<m_dof;i++) {
      for (int j=0;j<m_dof;j++) {
        m_val[i*m_dof + j] = 0.;
        for (int k=0;k<m_dof;k++) {
          m_val[i*m_dof + j] += a_scale.m_val[i*m_dof + k]*a_node.m_val[k*m_dof + j];
          if (m_val[i*m_dof + j] != m_val[i*m_dof + j]) {
            pout() << "StencilTensorValue:apply have a Nan: index:" <<i<<","<<j<<","<<k<<".  m_dof="<< m_dof<< std::endl;
            pout() << "\t\t scale is Nan:" << (a_scale.m_val[i*m_dof + k] != a_scale.m_val[i*m_dof + k]) << std::endl;
            pout() << "\t\t node is Nan:" << (a_node.m_val[k*m_dof + j] != a_node.m_val[k*m_dof + j]) << std::endl;
            MayDay::Error("PetscCompGrid::AddStencilToMatit->second.getVals(0,0) is a NaN");
          }
        }
      }
    }
  }
}
#include "NamespaceFooter.H"
