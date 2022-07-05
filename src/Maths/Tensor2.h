#pragma once

namespace almost
{
  template<typename T>
  class Tensor2
  {
  public:
    Tensor2()
    {
      xx = yy = xy = 0;
    }
    Tensor2(T xx, T xy, T yy)
    {
      this->xx = xx;
      this->xy = xy;
      this->yy = yy;
    }
    Tensor2(T diag)
    {
      xx = yy = diag;
      xy = 0;
    }
    glm::vec2 operator ()(const glm::vec2 &v) const
    {
      return glm::vec2(
        xx * v.x + xy * v.y,
        xy * v.x + yy * v.y);
    }
    T GetTrace()
    {
      return xx + yy;
    }

    Tensor2& operator +=(const T &v)
    {
      this->xx += v;
      this->yy += v;
      return *this;
    }

    Tensor2& operator +=(const Tensor2& t)
    {
      this->xx += t.xx;
      this->xy += t.xy;
      this->yy += t.yy;
      return *this;
    }

    Tensor2& operator *=(const T &v)
    {
      this->xx *= v;
      this->xy *= v;
      this->yy *= v;
      return *this;
    }

    Tensor2 GetInverse()
    {
      T det = xx * yy - xy * xy;
      return Tensor2(yy, -xy, xx) * (1.0f / det);
    }

    T xx, xy, yy;
  };

  template <typename T>
  Tensor2<T> TensorProduct(const glm::vec2& v)
  {
    return Tensor2<T>(
      v.x * v.x, //xx
      v.x * v.y, //xy
      v.y * v.y  //yy
    );
  }

  template <typename T>
  Tensor2<T> SymmetricTensorProduct(const glm::vec2 &v0, const glm::vec2 &v1) {
    return Tensor2<T>(
      T(2.0) * v0.x * v1.x, //xx
      v0.x * v1.y + v0.y * v1.x, //xy
      T(2.0) * v0.y * v1.y //yy
    );
  }


  template<typename T>
  inline const Tensor2<T> operator +(const Tensor2<T> &t0, const Tensor2<T> &t1)
  {
    return Tensor2<T>(
      t0.xx + t1.xx,
      t0.xy + t1.xy,
      t0.yy + t1.yy);
  }

  template<typename T>
  inline const Tensor2<T> operator +(const Tensor2<T> &t, const T v)
  {
    Tensor2<T> res = t;
    res += v;
    return res;
  }

  template<typename T>
  inline const Tensor2<T> operator -(const Tensor2<T> &t0, const Tensor2<T> &t1)
  {
    return Tensor2<T>(
      t0.xx - t1.xx,
      t0.xy - t1.xy,
      t0.yy - t1.yy);
  }

  template<typename T>
  inline const Tensor2<T> operator *(const Tensor2<T> &t, const T &s)
  {
    return Tensor2<T>(
      t.xx * s,
      t.xy * s,
      t.yy * s);
  }

  template<typename T>
  inline const Tensor2<T> operator *(const Tensor2<T> &t0, const Tensor2<T> &t1)
  {
    return Tensor2<T>(
      t0.xx * t1.xx + t0.xy * t1.xy,
      t0.xx * t1.xy + t0.xy * t1.yy,
      t0.xy * t1.xy + t0.yy * t1.yy);
  }

  template <class T>
  inline T DoubleConvolution(const Tensor2<T> &t0, const Tensor2<T> &t1)
  {
    return 
      t0.xx * t1.xx + 
      t0.yy * t1.yy + 
      T(2.0) * t0.xy * t1.xy;
  }

  typedef Tensor2<float>	Tensor2f;
  typedef Tensor2<double> Tensor2d;
}