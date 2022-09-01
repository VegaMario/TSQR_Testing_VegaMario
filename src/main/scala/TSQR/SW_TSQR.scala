package TSQR
import scala.collection.mutable
import scala.util.Random

object SW_TSQR {
  def main(args: Array[String]): Unit = {
    val ran = genRandomMatrix(6,3)
    ran.map(x=>println(x.map(_.round).toList))
    println()
//    val ran2 = genRandomMatrix(4,1)
//    ran2.map(x=>println(x.map(_.round).toList))
//    println()
    val result = Householder_QR_v2(ran)
    result.map(x=>println(x.map(_.round).toList))
    println()
  }
  def Householder_QR(matrix: Array[Array[Double]]): Array[Array[Double]] = {
    var matrix_cpy = matrix
    val m = matrix_cpy.length
    val n = matrix_cpy(0).length
    for(i <- 0 until n){
      val xk = getColumn(matrix_cpy,i,i,m)
      val vk = xk
      val xk_mag = magnitude(xk)
      vk(0) = xk(0) + sign(xk(0)) * xk_mag
      val vk_dot = dot_product(vk, vk)
      val t = -2/vk_dot
      for(j <- i until n){
        val xkj = getColumn(matrix_cpy,j,i,m)
        val xkj_mag = dot_product(xkj,vk)
        val tj = xkj_mag * t
        val tbs = axpy(tj,vk,xkj)
        matrix_cpy = update(matrix_cpy,tbs, j,i,m)
      }
      matrix_cpy.map(x=>println(x.map(_.round).toList))
      println()
    }
    matrix_cpy
  }

  def Householder_QR_v2(matrix: Array[Array[Double]]): Array[Array[Double]] = {
    val H_collection = mutable.ArrayBuffer[Array[Array[Double]]]()
    var matrix_cpy = matrix
    val m = matrix_cpy.length
    val n = matrix_cpy(0).length
    for(i <- 0 until n){
      val xk = getColumn(matrix_cpy,i,i,m)
      val vk = xk
      val xk_mag = magnitude(xk)
      vk(0) = xk(0) + sign(xk(0)) * xk_mag
      val vk_dot = dot_product(vk, vk)
      val vk_tp = tensorProduct(vk,vk)
      val I = genI(vk.length)
      val t = -2/vk_dot
      val nt = scalarMult(t,vk_tp)
      val H = matrixAdd(I,nt)
      H_collection += H
      val submatrix = matrixMultiply(H, getSubMatrix(matrix_cpy,i, n, i,m))
      matrix_cpy = replaceSM(matrix_cpy, submatrix,i,i)
      matrix_cpy.map(x=>println(x.map(_.round).toList))
      println()
    }
    val Q = Compute_Q(H_collection.toArray)
    for(i <- 0 until H_collection.length){
      println(s"m:${H_collection(i).length}")
      println(s"n:${H_collection(i)(0).length}")
      println()
    }
    println("Get original matrix from QR")
    val rebuild = matrixMultiply(Q, matrix_cpy)
    rebuild.map(x=>println(x.map(_.round).toList))
    println()
    matrix_cpy
  }
  def Compute_Q(Hs: Array[Array[Array[Double]]]): Array[Array[Double]] = {
    val fix = (for(i <- 0  until Hs.length)yield{
      val I = genI(Hs(0).length)
      val newa = replaceSM(I,Hs(i),i,i)
      newa
    }).map(_.map(_.toArray).toArray).toArray
    fix.map(x=>println(x.map(_.toList).toList)).toList
    val rev = fix.reverse
    var temp = rev(0)
    for(i <- 1 until rev.length){
      temp = matrixMultiply(rev(i), temp)
    }
    temp
  }
  def replaceSM(matrix: Array[Array[Double]], submatrix: Array[Array[Double]], col_ind: Int, start: Int): Array[Array[Double]] = {
    val matrix_cpy = matrix
    for(i <- col_ind until matrix_cpy(0).length){
      for(j <- start until matrix_cpy.length){
        matrix_cpy(j)(i) = submatrix(j-start)(i-col_ind)
      }
    }
    matrix_cpy
  }
  def scalarMult(scalar: Double, matrix: Array[Array[Double]]): Array[Array[Double]] = {
    matrix.map(_.map(_*scalar))
  }
  def matrixAdd(in_a: Array[Array[Double]], in_b: Array[Array[Double]]): Array[Array[Double]] = {
    val addition = (for(i <- 0 until in_a.length) yield{
      for(j <- 0 until in_a(0).length) yield{
        in_a(i)(j) + in_b(i)(j)
      }
    }).map(_.toArray).toArray
    addition
  }
  def matrixSub(in_a: Array[Array[Double]], in_b: Array[Array[Double]]): Array[Array[Double]] = {
    val subtraction = (for(i <- 0 until in_a.length) yield{
      for(j <- 0 until in_a(0).length) yield{
        in_a(i)(j) - in_b(i)(j)
      }
    }).map(_.toArray).toArray
    subtraction
  }
  def matrixMultiply(in_a: Array[Array[Double]], in_b: Array[Array[Double]]): Array[Array[Double]] = {
    val col_vecs = (for(i <- 0 until in_b(0).length) yield{
      for(j <- 0 until in_b.length)yield{
        in_b(j)(i)
      }
    }).map(_.toArray).toArray
    val multiplication = (for(i <- 0 until in_a.length)yield{
      for(j <-  0 until col_vecs.length) yield{
        dot_product(in_a(i), col_vecs(j))
      }
    }).map(_.toArray).toArray
    multiplication
  }
  def update(matrix: Array[Array[Double]], col: Array[Double], ind: Int, start: Int, stop: Int): Array[Array[Double]] = {
    val nm  = matrix
    for(i <- start until stop){
      nm(i)(ind) = col(i-start)
    }
    nm
  }
  def genI(n: Int): Array[Array[Double]] = {
    val I = (for (i <- 0 until n)yield{
      for(j <- 0 until n) yield{
        if(i == j)
          1.0
        else
          0.0
      }
    }).map(_.toArray).toArray
    I
  }
  def tensorProduct(in_a: Array[Double], in_b: Array[Double]): Array[Array[Double]] = {
    val newMatrix = (for(i <- 0 until in_a.length) yield{
      for(j <- 0 until in_b.length) yield{
        in_a(i) * in_b(j)
      }
    }).map(_.toArray).toArray
    newMatrix
  }
  def genRandomMatrix(m: Int, n: Int): Array[Array[Double]] = {
    val randmatrix = (for(i <- 0 until m) yield{
      for(j <- 0 until n) yield {
        val rand = new Random()
        val start = -100.0
        val end = 100.0
        (start + rand.nextLong((end - start).toLong + 1)) * 0.1
      }
    }).map(_.toArray).toArray
    randmatrix
  }
  def dot_product(in_a: Array[Double], in_b: Array[Double]): Double = {
    return in_a.zip(in_b).map(x=>x._1*x._2).reduce(_+_)
  }
  def sqrt(in_a: Double): Double = {
    return Math.sqrt(in_a)
  }
  def magnitude(in_a: Array[Double]): Double = {
    val dp = dot_product(in_a, in_a)
    sqrt(dp)
  }
  def sign(in_a: Double): Double = {
    if(in_a < -0.000000001)
      -1.0
    else
      1.0
  }
  def getColumn(matrix: Array[Array[Double]], ind: Int, start: Int, stop: Int): Array[Double] = {
    val column = (for(i <- start until stop) yield {
      matrix(i)(ind)
    }).toArray
    column
  }
  def getSubMatrix(matrix: Array[Array[Double]], ind_start: Int, ind_stop: Int, start: Int, stop: Int): Array[Array[Double]] = {
    val submatrix = (for(i <- start until stop) yield{
      for(j <- ind_start until ind_stop) yield{
        matrix(i)(j)
      }
    }).map(_.toArray).toArray
    submatrix
  }
  def axpy(a: Double, x: Array[Double], y: Array[Double]): Array[Double] = {
    val ax = x.map(x=>x*a)
    val axy = ax.zip(y).map(x=>x._1+x._2)
    axy
  }
}
