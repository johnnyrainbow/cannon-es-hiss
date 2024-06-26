import { Equation } from '../equations/Equation'
import { Vec3 } from '../math/Vec3'
import type { Body } from '../objects/Body'

/**
 * Contact/non-penetration constraint equation
 */
export class ContactEquation extends Equation {
  /**
   * "bounciness": u1 = -e*u0
   */
  restitution: number
  /**
   * World-oriented vector that goes from the center of bi to the contact point.
   */
  ri: Vec3
  /**
   * World-oriented vector that starts in body j position and goes to the contact point.
   */
  rj: Vec3
  /**
   * Contact normal, pointing out of body i.
   */
  ni: Vec3

  constructor(bodyA: Body, bodyB: Body, maxForce = 1e6) {
    super(bodyA, bodyB, 0, maxForce)

    this.restitution = 0.0
    this.ri = new Vec3()
    this.rj = new Vec3()
    this.ni = new Vec3()
  }

  computeB(h: number): number {
    // console.log('computing b CONTACT')

    const a = this.a
    const b = this.b
    const bi = this.bi
    const bj = this.bj
    const ri = this.ri
    const rj = this.rj
    const rixn = ContactEquation_computeB_temp1
    const rjxn = ContactEquation_computeB_temp2
    let vi = bi.velocity
    const wi = bi.angularVelocity
    const fi = bi.force
    const taui = bi.torque
    let vj = bj.velocity
    const wj = bj.angularVelocity
    wj.x = 0
    wj.y = 0
    wj.z = 0
    wi.x = 0
    wi.y = 0
    wi.z = 0

    // vi.x = 0
    // vi.y = 0
    // vi.z = 0

    // if (bi.customIndependentVelocity.x !== 0 || bi.customIndependentVelocity.z !== 0) {
    //   vi = bi.customIndependentVelocity
    // }
    // if (bj.customIndependentVelocity.x !== 0 || bj.customIndependentVelocity.z !== 0) {
    //   vj = bj.customIndependentVelocity
    // }

    const fj = bj.force
    const tauj = bj.torque
    const penetrationVec = ContactEquation_computeB_temp3
    const GA = this.jacobianElementA
    const GB = this.jacobianElementB
    const n = this.ni

    // Caluclate cross products
    ri.cross(n, rixn)
    rj.cross(n, rjxn)

    // g = xj+rj -(xi+ri)
    // G = [ -ni  -rixn  ni  rjxn ]
    n.negate(GA.spatial)
    rixn.negate(GA.rotational)
    GB.spatial.copy(n)
    GB.rotational.copy(rjxn)

    // Calculate the penetration vector
    const bjPC = new Vec3().copy(bj.position)
    const biPC = new Vec3().copy(bi.position)
    // bjPC.y = 0
    // biPC.y = 0
    // console.log('used pen pre', penetrationVec, rj, ri)
    penetrationVec.copy(bjPC)
    penetrationVec.vadd(rj, penetrationVec)
    penetrationVec.vsub(biPC, penetrationVec)
    penetrationVec.vsub(ri, penetrationVec)
    // penetrationVec.x = 0
    // penetrationVec.y = 0
    // penetrationVec.z = 0
    let g = n.dot(penetrationVec)

    // Compute iteration
    const ePlusOne = this.restitution + 1
    const GW = ePlusOne * vj.dot(n) - ePlusOne * vi.dot(n) + wj.dot(rjxn) - wi.dot(rixn)
    const GiMf = this.computeGiMf()

    //TODO variants are -g, GW, giMF
    // console.log('used n', n)
    // console.log('used pen', penetrationVec)
    // console.log('used pos', bjPC, biPC)
    // // g = -2
    // console.log('used ni', -g, a, GW, b, h, GiMf)
    const B = -g * a - GW * b - h * GiMf

    return B
  }

  /**
   * Get the current relative velocity in the contact point.
   */
  getImpactVelocityAlongNormal(): number {
    const vi = ContactEquation_getImpactVelocityAlongNormal_vi
    const vj = ContactEquation_getImpactVelocityAlongNormal_vj
    const xi = ContactEquation_getImpactVelocityAlongNormal_xi
    const xj = ContactEquation_getImpactVelocityAlongNormal_xj
    const relVel = ContactEquation_getImpactVelocityAlongNormal_relVel

    this.bi.position.vadd(this.ri, xi)
    this.bj.position.vadd(this.rj, xj)

    this.bi.getVelocityAtWorldPoint(xi, vi)
    this.bj.getVelocityAtWorldPoint(xj, vj)

    vi.vsub(vj, relVel)
    console.log('RELVEL', relVel)
    return this.ni.dot(relVel)
  }
}

const ContactEquation_computeB_temp1 = new Vec3() // Temp vectors
const ContactEquation_computeB_temp2 = new Vec3()
const ContactEquation_computeB_temp3 = new Vec3()

const ContactEquation_getImpactVelocityAlongNormal_vi = new Vec3()
const ContactEquation_getImpactVelocityAlongNormal_vj = new Vec3()
const ContactEquation_getImpactVelocityAlongNormal_xi = new Vec3()
const ContactEquation_getImpactVelocityAlongNormal_xj = new Vec3()
const ContactEquation_getImpactVelocityAlongNormal_relVel = new Vec3()
