#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif
#ifdef CH_USE_PETSC

#include "PetscAMRSolver.H"
#include "JFNKSolver.H"
#include "NamespaceHeader.H"

PetscAMRSolver::PetscAMRSolver() : m_op_mfree(0),m_mfree_homogeneous(true)
{
}
#undef __FUNCT__
#define __FUNCT__ "apply_mfree"
PetscErrorCode PetscAMRSolver::apply_mfree(Mat A, Vec x, Vec f)
{
  CH_TIME("PetscAMRSolver::apply_mfree");
  //Whenever I write any PETSc code, I look forward to pulling classes back from the void.
  PetscFunctionBegin;
  void *ctx;
  MatShellGetContext(A, &ctx);
  PetscAMRSolver *tthis = (PetscAMRSolver*)ctx;
  tthis->m_petscCompMat.putPetscInChombo(x, tthis->m_phi_mfree);
  tthis->m_op_mfree->applyOp(tthis->m_Lphi_mfree,tthis->m_phi_mfree,tthis->m_mfree_homogeneous);
  tthis->m_petscCompMat.putChomboInPetsc(tthis->m_Lphi_mfree,f);
  PetscFunctionReturn(0);
}
#undef __FUNCT__
#define __FUNCT__ "solve_mfree"
PetscErrorCode
PetscAMRSolver::solve_mfree( Vector<LevelData<FArrayBox>*>& a_phi, 
			     const Vector<LevelData<FArrayBox>*>& a_rhs, 
			     JFNKOp *a_op )
{
  CH_TIME("PetscSolver::solve_mfree");
#ifdef CH_MPI
  MPI_Comm wcomm = Chombo_MPI::comm;
#else
  MPI_Comm wcomm = PETSC_COMM_SELF;
#endif
  PetscErrorCode ierr;
  Vec  x, b;      /* approx solution, RHS */
  Mat  A;         /* linear system matrix */
  KSP  ksp;       /* linear solver context */
  PetscFunctionBeginUser;

  ierr = m_petscCompMat.createMatrix(); CHKERRQ(ierr); 
  A = m_petscCompMat.getMatrix();

  //create an operator matrix shell with same dimensions as m_mat
  PetscInt m, n, M, N;
  ierr = MatGetSize(A, &M, &N);CHKERRQ(ierr); CH_assert(M == N);
  ierr = MatGetLocalSize(A, &m, &n);CHKERRQ(ierr);
  Mat L; 
  ierr = MatCreateShell(wcomm,m,n,N,N,(void *)this,&L);CHKERRQ(ierr);
  ierr = MatShellSetOperation(L,MATOP_MULT,(void(*)(void))apply_mfree);
  m_op_mfree = a_op;
  //allocate space for a vector and a matrix-vector product in Chombo-land
  a_op->create( m_phi_mfree , a_phi);
  a_op->create( m_Lphi_mfree , a_rhs);

  ierr = MatGetVecs(A,&x,&b); CHKERRQ(ierr);
  ierr = m_petscCompMat.putChomboInPetsc(a_rhs,b); CHKERRQ(ierr);
  ierr = KSPCreate(wcomm, &ksp); CHKERRQ(ierr);
  ierr = KSPSetOperators(ksp, L, A, SAME_NONZERO_PATTERN); CHKERRQ(ierr);
  ierr = KSPSetFromOptions(ksp); CHKERRQ(ierr);
  ierr = KSPSolve(ksp, b, x); CHKERRQ(ierr);
  
  ierr = m_petscCompMat.putPetscInChombo(x,a_phi); CHKERRQ(ierr);
  
  ierr = VecDestroy(&x); CHKERRQ(ierr);
  ierr = VecDestroy(&b); CHKERRQ(ierr);
  ierr = KSPDestroy(&ksp); CHKERRQ(ierr);
  ierr = MatDestroy(&L); CHKERRQ(ierr);

  //clean up 
  a_op->clear(m_phi_mfree);
  a_op->clear(m_Lphi_mfree);
  PetscFunctionReturn(0); 
}
#undef __FUNCT__
#define __FUNCT__ "solve"
PetscErrorCode
PetscAMRSolver::solve( Vector<LevelData<FArrayBox>*>& a_phi, 
		       const Vector<LevelData<FArrayBox>*>& a_rhs )
{
  CH_TIME("PetscSolver::solve_mfree");
#ifdef CH_MPI
  MPI_Comm wcomm = Chombo_MPI::comm;
#else
  MPI_Comm wcomm = PETSC_COMM_SELF;
#endif
  PetscErrorCode ierr;
  Vec  x, b;      /* approx solution, RHS */
  Mat  A;         /* linear system matrix */
  KSP  ksp;       /* linear solver context */
  PetscFunctionBeginUser;

  ierr = m_petscCompMat.createMatrix(); CHKERRQ(ierr); 
  A = m_petscCompMat.getMatrix();
  
  //create an operator matrix shell with same dimensions as m_mat
  ierr = MatGetVecs(A,&x,&b); CHKERRQ(ierr);
  ierr = m_petscCompMat.putChomboInPetsc(a_rhs,b); CHKERRQ(ierr);
  ierr = KSPCreate(wcomm, &ksp); CHKERRQ(ierr);
  ierr = KSPSetOperators(ksp, A, A, SAME_NONZERO_PATTERN); CHKERRQ(ierr);
  ierr = KSPSetFromOptions(ksp); CHKERRQ(ierr);
  ierr = KSPSolve(ksp, b, x); CHKERRQ(ierr);
  
  ierr = m_petscCompMat.putPetscInChombo(x,a_phi); CHKERRQ(ierr);
  
  ierr = VecDestroy(&x); CHKERRQ(ierr);
  ierr = VecDestroy(&b); CHKERRQ(ierr);
  ierr = KSPDestroy(&ksp); CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

#include "NamespaceFooter.H"

#endif
