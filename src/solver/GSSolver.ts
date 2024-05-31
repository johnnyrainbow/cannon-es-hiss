import { Solver } from '../solver/Solver'
import type { World } from '../world/World'

/**
 * Constraint equation Gauss-Seidel solver.
 * @todo The spook parameters should be specified for each constraint, not globally.
 * @see https://www8.cs.umu.se/kurser/5DV058/VT09/lectures/spooknotes.pdf
 */
export class GSSolver extends Solver {
  /**
   * The number of solver iterations determines quality of the constraints in the world.
   * The more iterations, the more correct simulation. More iterations need more computations though. If you have a large gravity force in your world, you will need more iterations.
   */
  iterations: number

  /**
   * When tolerance is reached, the system is assumed to be converged.
   */
  tolerance: number

  /**
   * @todo remove useless constructor
   */
  constructor() {
    super()

    this.iterations = 10
    this.tolerance = 1e-7
  }

  /**
   * Solve
   * @return number of iterations performed
   */
  solve(dt: number, world: World): number {
    let iter = 0
    const maxIter = this.iterations
    const tolSquared = this.tolerance * this.tolerance
    const equations = this.equations
    const Neq = equations.length
    const bodies = world.bodies
    const Nbodies = bodies.length
    const h = dt
    let q
    let B
    let invC
    let deltalambda
    let deltalambdaTot
    let GWlambda
    let lambdaj

    // Update solve mass
    if (Neq !== 0) {
      for (let i = 0; i !== Nbodies; i++) {
        bodies[i].updateSolveMassProperties()
      }
    }

    // Things that do not change during iteration can be computed once
    const invCs = GSSolver_solve_invCs

    const Bs = GSSolver_solve_Bs
    const lambda = GSSolver_solve_lambda
    invCs.length = Neq
    Bs.length = Neq
    lambda.length = Neq
    for (let i = 0; i !== Neq; i++) {
      const c = equations[i] as any
      lambda[i] = 0.0
      Bs[i] = c.computeB(h)
      invCs[i] = 1.0 / c.computeC()
    }

    if (Neq !== 0) {
      // Reset vlambda
      for (let i = 0; i !== Nbodies; i++) {
        const b = bodies[i]
        const vlambda = b.vlambda
        const wlambda = b.wlambda
        vlambda.set(0, 0, 0)
        wlambda.set(0, 0, 0)
      }

      // Iterate over equations
      for (iter = 0; iter !== maxIter; iter++) {
        // Accumulate the total error for each iteration.
        deltalambdaTot = 0.0

        for (let j = 0; j !== Neq; j++) {
          const c = equations[j]

          // Compute iteration
          B = Bs[j]
          invC = invCs[j]
          lambdaj = lambda[j]
          GWlambda = c.computeGWlambda()
          console.log('computed lambda', GWlambda)
          deltalambda = invC * (B - GWlambda - c.eps * lambdaj)

          // Clamp if we are not within the min/max interval
          if (lambdaj + deltalambda < c.minForce) {
            deltalambda = c.minForce - lambdaj
          } else if (lambdaj + deltalambda > c.maxForce) {
            deltalambda = c.maxForce - lambdaj
          }
          lambda[j] += deltalambda

          deltalambdaTot += deltalambda > 0.0 ? deltalambda : -deltalambda // abs(deltalambda)
          if (deltalambda > 50) {
            // console.log(invC, B, GWlambda, c.eps, lambdaj)
          }
          c.addToWlambda(deltalambda)
        }

        // If the total error is small enough - stop iterate
        if (deltalambdaTot * deltalambdaTot < tolSquared) {
          break
        }
      }

      // Add result to velocity
      for (let i = 0; i !== Nbodies; i++) {
        const b = bodies[i]

        //GAB ITS HERE!!!!

        //it seems the velocity will be mutated by the collision shit.
        //so if we hit a wall, the velocity will be set to 0 here.
        //therefore, we must apply the mutation to the customIndependentVelocity also.
        //so that collisions will work.
        const v = b.velocity
        const v2 = b.customIndependentVelocity
        const w = b.angularVelocity
        b.vlambda.vmul(b.linearFactor, b.vlambda)

        v.vadd(b.vlambda, v)
        const posv2x = Math.abs(b.vlambda.x)
        const posv2z = Math.abs(b.vlambda.z)
        // b.vlambda.x = parseFloat(b.vlambda.x.toFixed(1))
        // b.vlambda.y = parseFloat(b.vlambda.x.toFixed(1))
        // b.vlambda.z = parseFloat(b.vlambda.x.toFixed(1))
        // v2.vadd(b.vlambda, v2)
        // v2.vadd(b.vlambda, v2)
        let modd = false
        if (Math.abs(b.vlambda.x) > 400) {
          console.log('PAPA2')
          console.log('DID you knowww x', b.position, b.customIndependentVelocity, b.vlambda.x)
          if (b.vlambda.x < 0) {
            b.vlambda.x = -2000
          } else {
            b.vlambda.x = 2000
          }
          modd = true
        } else {
          b.vlambda.x = 0
        }
        //test2
        if (Math.abs(b.vlambda.z) > 400) {
          console.log('DID you knowww z', b.position, b.customIndependentVelocity, b.vlambda.z)
          if (b.vlambda.z < 0) {
            b.vlambda.z = -2000
          } else {
            b.vlambda.z = 2000
          }
          modd = true
        } else {
          b.vlambda.z = 0
        }

        // modd = true;
        if (modd) {
          v2.vadd(b.vlambda, v2)
          // if (Math.abs(b.vlambda.x) > 100 || Math.abs(b.vlambda.z) > 100) {

          console.log('DID YOU KNOW YOUR vlambdaoo7?', b.vlambda)
          // }
          //F
        }

        // if (Math.abs(v2.x) > 30 || Math.abs(v2.z) > 30) {
        //   // console.log("SET CUSTOM VELOOCH ", v2)
        //   v2.vadd(b.vlambda, v2)
        // }
        b.wlambda.vmul(b.angularFactor, b.wlambda)
        w.vadd(b.wlambda, w)
      }

      // Set the `.multiplier` property of each equation
      let l = equations.length
      const invDt = 1 / h
      while (l--) {
        equations[l].multiplier = lambda[l] * invDt
      }
    }

    return iter
  }
}

// Just temporary number holders that we want to reuse each iteration.
const GSSolver_solve_lambda: number[] = []
const GSSolver_solve_invCs: number[] = []
const GSSolver_solve_Bs: number[] = []
