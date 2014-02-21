/*
 * Copyright 2014 MovingBlocks
 *
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 *
 *  http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 */

package org.terasology.math.geom;

/**
 * A double precision floating point 3 by 3 matrix.
 * @author Martin Steiger
 */
public class ArrayBasedMatrix3d extends BaseMatrix3d {

    private double[] arr = new double[9];

    /**
     * Constructs and initializes a Matrix3d from the specified nine values.
     * @param m00 the [0][0] element
     * @param m01 the [0][1] element
     * @param m02 the [0][2] element
     * @param m10 the [1][0] element
     * @param m11 the [1][1] element
     * @param m12 the [1][2] element
     * @param m20 the [2][0] element
     * @param m21 the [2][1] element
     * @param m22 the [2][2] element
     */
    public ArrayBasedMatrix3d(double m00, double m01, double m02,
            double m10, double m11, double m12,
            double m20, double m21, double m22) {
        this.set(0, 0, m00);
        this.set(0, 1, m01);
        this.set(0, 2, m02);

        this.set(1, 0, m10);
        this.set(1, 1, m11);
        this.set(1, 2, m12);

        this.set(2, 0, m20);
        this.set(2, 1, m21);
        this.set(2, 2, m22);

    }

    /**
     * Constructs and initializes a Matrix3d from the specified nine-
     * element array.
     * @param v the array of length 9 containing in order
     */
    public ArrayBasedMatrix3d(double[] v) {
        this(v[0], v[1], v[2],
             v[3], v[4], v[5],
             v[6], v[7], v[8]);
    }

    /**
      *  Constructs a new matrix with the same values as the
      *  Matrix3d parameter.
      *  @param m1  the source matrix
      */
    public ArrayBasedMatrix3d(ArrayBasedMatrix3d m1) {
        this(
            m1.getM00(),
            m1.getM01(),
            m1.getM02(),
    
            m1.getM10(),
            m1.getM11(),
            m1.getM12(),
    
            m1.getM20(),
            m1.getM21(),
            m1.getM22());
    }

    /**
     * Constructs and initializes a Matrix3d to all zeros.
     */
    public ArrayBasedMatrix3d() {
        this.set(0, 0, 0.0);
        this.set(0, 1, 0.0);
        this.set(0, 2, 0.0);

        this.set(1, 0, 0.0);
        this.set(1, 1, 0.0);
        this.set(1, 2, 0.0);

        this.set(2, 0, 0.0);
        this.set(2, 1, 0.0);
        this.set(2, 2, 0.0);
    }

    /**
     * Sets this Matrix3d to identity.
     */
    public final void setIdentity() {
        this.set(0, 0, 1.0);
        this.set(0, 1, 0.0);
        this.set(0, 2, 0.0);

        this.set(1, 0, 0.0);
        this.set(1, 1, 1.0);
        this.set(1, 2, 0.0);

        this.set(2, 0, 0.0);
        this.set(2, 1, 0.0);
        this.set(2, 2, 1.0);
    }

    /**
     * Retrieves the value at the specified row and column of the specified
     * matrix.
     * @param row the row number to be retrieved (zero indexed)
     * @param column the column number to be retrieved (zero indexed)
     * @return the value at the indexed element.
     */
    @Override
    public final double get(int row, int column) {
        if (row < 0 || row >= 3) {
            throw new ArrayIndexOutOfBoundsException("col not in [0..2]");
        }

        if (column < 0 || column >= 3) {
            throw new ArrayIndexOutOfBoundsException("col not in [0..2]");
        }

        return arr[row * 3 + column];
    }

    /**
     * Sets the specified element of this matrix3f to the value provided.
     * @param row the row number to be modified (zero indexed)
     * @param column the column number to be modified (zero indexed)
     * @param value the new value
     */
    public final void set(int row, int column, double value) {
        if (row < 0 || row >= 3) {
            throw new ArrayIndexOutOfBoundsException("col not in [0..2]");
        }

        if (column < 0 || column >= 3) {
            throw new ArrayIndexOutOfBoundsException("col not in [0..2]");
        }

        arr[row * 3 + column] = value;
    }

    /**
     * Sets the specified row of this matrix3d to the 4 values provided.
     * @param row the row number to be modified (zero indexed)
     * @param x the first column element
     * @param y the second column element
     * @param z the third column element
     */
    public final void setRow(int row, double x, double y, double z) {
        switch (row) {
            case 0:
                this.set(0, 0, x);
                this.set(0, 1, y);
                this.set(0, 2, z);
                break;

            case 1:
                this.set(1, 0, x);
                this.set(1, 1, y);
                this.set(1, 2, z);
                break;

            case 2:
                this.set(2, 0, x);
                this.set(2, 1, y);
                this.set(2, 2, z);
                break;

            default:
                throw new ArrayIndexOutOfBoundsException("row not in [0..2]");
        }
    }

    /**
     * Sets the specified row of this matrix3d to the Vector provided.
     * @param row the row number to be modified (zero indexed)
     * @param v the replacement row
     */
    public final void setRow(int row, Vector3d v) {
        switch (row) {
            case 0:
                this.set(0, 0, v.getX());
                this.set(0, 1, v.getY());
                this.set(0, 2, v.getZ());
                break;

            case 1:
                this.set(1, 0, v.getX());
                this.set(1, 1, v.getY());
                this.set(1, 2, v.getZ());
                break;

            case 2:
                this.set(2, 0, v.getX());
                this.set(2, 1, v.getY());
                this.set(2, 2, v.getZ());
                break;

            default:
                throw new ArrayIndexOutOfBoundsException("row not in [0..2]");
        }
    }

    /**
     * Sets the specified row of this matrix3d to the three values provided.
     * @param row the row number to be modified (zero indexed)
     * @param v the replacement row
     */
    public final void setRow(int row, double[] v) {
        switch (row) {
            case 0:
                this.set(0, 0, v[0]);
                this.set(0, 1, v[1]);
                this.set(0, 2, v[2]);
                break;

            case 1:
                this.set(1, 0, v[0]);
                this.set(1, 1, v[1]);
                this.set(1, 2, v[2]);
                break;

            case 2:
                this.set(2, 0, v[0]);
                this.set(2, 1, v[1]);
                this.set(2, 2, v[2]);
                break;

            default:
                throw new ArrayIndexOutOfBoundsException("row not in [0..2]");
        }
    }

    /**
     * Sets the specified column of this matrix3d to the three values provided.
     * @param column the column number to be modified (zero indexed)
     * @param x the first row element
     * @param y the second row element
     * @param z the third row element
     */
    public final void setColumn(int column, double x, double y, double z) {
        switch (column) {
            case 0:
                this.set(0, 0, x);
                this.set(1, 0, y);
                this.set(2, 0, z);
                break;

            case 1:
                this.set(0, 1, x);
                this.set(1, 1, y);
                this.set(2, 1, z);
                break;

            case 2:
                this.set(0, 2, x);
                this.set(1, 2, y);
                this.set(2, 2, z);
                break;

            default:
                throw new ArrayIndexOutOfBoundsException("col not in [0..2]");
        }
    }

    /**
     * Sets the specified column of this matrix3d to the vector provided.
     * @param column the column number to be modified (zero indexed)
     * @param v the replacement column
     */
    public final void setColumn(int column, Vector3d v) {
        switch (column) {
            case 0:
                this.set(0, 0, v.getX());
                this.set(1, 0, v.getY());
                this.set(2, 0, v.getZ());
                break;

            case 1:
                this.set(0, 1, v.getX());
                this.set(1, 1, v.getY());
                this.set(2, 1, v.getZ());
                break;

            case 2:
                this.set(0, 2, v.getX());
                this.set(1, 2, v.getY());
                this.set(2, 2, v.getZ());
                break;

            default:
                throw new ArrayIndexOutOfBoundsException("col not in [0..2]");
        }
    }

    /**
     * Sets the specified column of this matrix3d to the three values provided.
     * @param column the column number to be modified (zero indexed)
     * @param v the replacement column
     */
    public final void setColumn(int column, double[] v) {
        switch (column) {
            case 0:
                this.set(0, 0, v[0]);
                this.set(1, 0, v[1]);
                this.set(2, 0, v[2]);
                break;

            case 1:
                this.set(0, 1, v[0]);
                this.set(1, 1, v[1]);
                this.set(2, 1, v[2]);
                break;

            case 2:
                this.set(0, 2, v[0]);
                this.set(1, 2, v[1]);
                this.set(2, 2, v[2]);
                break;

            default:
                throw new ArrayIndexOutOfBoundsException("col not in [0..2]");
        }
    }

    /**
     * Adds a scalar to each component of this matrix.
     * @param scalar  the scalar adder
     */
    public final void add(double scalar) {
        set(0, 0, get(0, 0) + scalar);
        set(0, 1, get(0, 1) + scalar);
        set(0, 2, get(0, 2) + scalar);

        set(1, 0, get(1, 0) + scalar);
        set(1, 1, get(1, 1) + scalar);
        set(1, 2, get(1, 2) + scalar);

        set(2, 0, get(2, 0) + scalar);
        set(2, 1, get(2, 1) + scalar);
        set(2, 2, get(2, 2) + scalar);

    }

    /**
     * Sets the value of this matrix to the sum of itself and matrix m1.
     * @param m1 the other matrix
     */
    public final void add(ArrayBasedMatrix3d m1) {
        this.set(0, 0, get(0, 0) + m1.get(0, 0));
        this.set(0, 1, get(0, 1) + m1.get(0, 1));
        this.set(0, 2, get(0, 2) + m1.get(0, 2));

        this.set(1, 0, get(1, 0) + m1.get(1, 0));
        this.set(1, 1, get(1, 1) + m1.get(1, 1));
        this.set(1, 2, get(1, 2) + m1.get(1, 2));

        this.set(2, 0, get(2, 0) + m1.get(2, 0));
        this.set(2, 1, get(2, 1) + m1.get(2, 1));
        this.set(2, 2, get(2, 2) + m1.get(2, 2));
    }

    /**
     * Sets the value of this matrix to the matrix difference of itself and
     * matrix m1 (this = this - m1).
     * @param m1 the other matrix
     */
    public final void sub(ArrayBasedMatrix3d m1) {
        this.set(0, 0, get(0, 0) - m1.get(0, 0));
        this.set(0, 1, get(0, 1) - m1.get(0, 1));
        this.set(0, 2, get(0, 2) - m1.get(0, 2));

        this.set(1, 0, get(1, 0) - m1.get(1, 0));
        this.set(1, 1, get(1, 1) - m1.get(1, 1));
        this.set(1, 2, get(1, 2) - m1.get(1, 2));

        this.set(2, 0, get(2, 0) - m1.get(2, 0));
        this.set(2, 1, get(2, 1) - m1.get(2, 1));
        this.set(2, 2, get(2, 2) - m1.get(2, 2));
    }

    /**
     * Sets the value of this matrix to its transpose.
     */
    public final void transpose() {
        double temp;

        temp = this.get(1, 0);
        this.set(1, 0, this.get(0, 1));
        this.set(0, 1, temp);

        temp = this.get(2, 0);
        this.set(2, 0, this.get(0, 2));
        this.set(0, 2, temp);

        temp = this.get(2, 1);
        this.set(2, 1, this.get(1, 2));
        this.set(1, 2, temp);
    }

    /**
     * Sets the value of this matrix to the transpose of the argument matrix.
     * @param m1 the matrix to be transposed
     */
    public final void transpose(ArrayBasedMatrix3d m1) {
        this.set(0, 0, m1.get(0, 0));
        this.set(0, 1, m1.get(1, 0));
        this.set(0, 2, m1.get(2, 0));

        this.set(1, 0, m1.get(0, 1));
        this.set(1, 1, m1.get(1, 1));
        this.set(1, 2, m1.get(2, 1));

        this.set(2, 0, m1.get(0, 2));
        this.set(2, 1, m1.get(1, 2));
        this.set(2, 2, m1.get(2, 2));
    }

    /**
     * Sets the value of this matrix to the matrix conversion of the
     * double precision quaternion argument.
     * @param q1 the quaternion to be converted
     */
    public final void set(Quat4d q1) {
        this.set(0, 0, (1.0 - 2.0 * q1.getY() * q1.getY() - 2.0 * q1.getZ() * q1.getZ()));
        this.set(1, 0, (2.0 * (q1.getX() * q1.getY() + q1.getW() * q1.getZ())));
        this.set(2, 0, (2.0 * (q1.getX() * q1.getZ() - q1.getW() * q1.getY())));

        this.set(0, 1, (2.0 * (q1.getX() * q1.getY() - q1.getW() * q1.getZ())));
        this.set(1, 1, (1.0 - 2.0 * q1.getX() * q1.getX() - 2.0 * q1.getZ() * q1.getZ()));
        this.set(2, 1, (2.0 * (q1.getY() * q1.getZ() + q1.getW() * q1.getX())));

        this.set(0, 2, (2.0 * (q1.getX() * q1.getZ() + q1.getW() * q1.getY())));
        this.set(1, 2, (2.0 * (q1.getY() * q1.getZ() - q1.getW() * q1.getX())));
        this.set(2, 2, (1.0 - 2.0 * q1.getX() * q1.getX() - 2.0 * q1.getY() * q1.getY()));
    }

    /**
     * Sets the value of this matrix to the value of the Matrix3d
     * argument.
     * @param m1 the source matrix3d
     */
    public final void set(ArrayBasedMatrix3d m1) {
        this.set(0, 0, m1.get(0, 0));
        this.set(0, 1, m1.get(0, 1));
        this.set(0, 2, m1.get(0, 2));

        this.set(1, 0, m1.get(1, 0));
        this.set(1, 1, m1.get(1, 1));
        this.set(1, 2, m1.get(1, 2));

        this.set(2, 0, m1.get(2, 0));
        this.set(2, 1, m1.get(2, 1));
        this.set(2, 2, m1.get(2, 2));
    }

    /**
     *  Sets the values in this Matrix3d equal to the row-major
     *  array parameter (ie, the first three elements of the
     *  array will be copied into the first row of this matrix, etc.).
     *  @param m  the double precision array of length 9
     */
    public final void set(double[] m) {
        set(0, 0, m[0]);
        set(0, 1, m[1]);
        set(0, 2, m[2]);

        set(1, 0, m[3]);
        set(1, 1, m[4]);
        set(1, 2, m[5]);

        set(2, 0, m[6]);
        set(2, 1, m[7]);
        set(2, 2, m[8]);

    }

    /**
     * Sets the value of this matrix to a scale matrix with
     * the passed scale amount.
     * @param scale the scale factor for the matrix
     */
    public final void set(double scale) {
        this.set(0, 0, scale);
        this.set(0, 1, 0.0);
        this.set(0, 2, 0.0);

        this.set(1, 0, 0.0);
        this.set(1, 1, scale);
        this.set(1, 2, 0.0);

        this.set(2, 0, 0.0);
        this.set(2, 1, 0.0);
        this.set(2, 2, scale);
    }

    /**
     * Sets the value of this matrix to a counter clockwise rotation
     * about the x axis.
     * @param angle the angle to rotate about the X axis in radians
     */
    public final void setRotX(double angle) {
        double sinAngle;
        double cosAngle;

        sinAngle = Math.sin(angle);
        cosAngle = Math.cos(angle);

        this.set(0, 0, 1.0);
        this.set(0, 1, 0.0);
        this.set(0, 2, 0.0);

        this.set(1, 0, 0.0);
        this.set(1, 1, cosAngle);
        this.set(1, 2, -sinAngle);

        this.set(2, 0, 0.0);
        this.set(2, 1, sinAngle);
        this.set(2, 2, cosAngle);
    }

    /**
     * Sets the value of this matrix to a counter clockwise rotation
     * about the y axis.
     * @param angle the angle to rotate about the Y axis in radians
     */
    public final void setRotY(double angle) {
        double sinAngle;
        double cosAngle;

        sinAngle = Math.sin(angle);
        cosAngle = Math.cos(angle);

        this.set(0, 0, cosAngle);
        this.set(0, 1, 0.0);
        this.set(0, 2, sinAngle);

        this.set(1, 0, 0.0);
        this.set(1, 1, 1.0);
        this.set(1, 2, 0.0);

        this.set(2, 0, -sinAngle);
        this.set(2, 1, 0.0);
        this.set(2, 2, cosAngle);
    }

    /**
     * Sets the value of this matrix to a counter clockwise rotation
     * about the z axis.
     * @param angle the angle to rotate about the Z axis in radians
     */
    public final void rotZ(double angle) {
        double sinAngle;
        double cosAngle;

        sinAngle = Math.sin(angle);
        cosAngle = Math.cos(angle);

        this.set(0, 0, cosAngle);
        this.set(0, 1, -sinAngle);
        this.set(0, 2, 0.0);

        this.set(1, 0, sinAngle);
        this.set(1, 1, cosAngle);
        this.set(1, 2, 0.0);

        this.set(2, 0, 0.0);
        this.set(2, 1, 0.0);
        this.set(2, 2, 1.0);
    }

    /**
      * Multiplies each element of this matrix by a scalar.
      * @param scalar  The scalar multiplier.
      */
    public final void mul(double scalar) {
        set(0, 0, get(0, 0) * scalar);
        set(0, 1, get(0, 1) * scalar);
        set(0, 2, get(0, 2) * scalar);

        set(1, 0, get(1, 0) * scalar);
        set(1, 1, get(1, 1) * scalar);
        set(1, 2, get(1, 2) * scalar);

        set(2, 0, get(2, 0) * scalar);
        set(2, 1, get(2, 1) * scalar);
        set(2, 2, get(2, 2) * scalar);

    }

    /**
      * Sets the value of this matrix to the result of multiplying itself
      * with matrix m1.
      * @param m1 the other matrix
      */
    public final void mul(ArrayBasedMatrix3d m1) {
        double lm00;
        double lm01;
        double lm02;
        double lm10;
        double lm11;
        double lm12;
        double lm20;
        double lm21;
        double lm22;

        lm00 = this.get(0, 0) * m1.get(0, 0) + this.get(0, 1) * m1.get(1, 0) + this.get(0, 2) * m1.get(2, 0);
        lm01 = this.get(0, 0) * m1.get(0, 1) + this.get(0, 1) * m1.get(1, 1) + this.get(0, 2) * m1.get(2, 1);
        lm02 = this.get(0, 0) * m1.get(0, 2) + this.get(0, 1) * m1.get(1, 2) + this.get(0, 2) * m1.get(2, 2);

        lm10 = this.get(1, 0) * m1.get(0, 0) + this.get(1, 1) * m1.get(1, 0) + this.get(1, 2) * m1.get(2, 0);
        lm11 = this.get(1, 0) * m1.get(0, 1) + this.get(1, 1) * m1.get(1, 1) + this.get(1, 2) * m1.get(2, 1);
        lm12 = this.get(1, 0) * m1.get(0, 2) + this.get(1, 1) * m1.get(1, 2) + this.get(1, 2) * m1.get(2, 2);

        lm20 = this.get(2, 0) * m1.get(0, 0) + this.get(2, 1) * m1.get(1, 0) + this.get(2, 2) * m1.get(2, 0);
        lm21 = this.get(2, 0) * m1.get(0, 1) + this.get(2, 1) * m1.get(1, 1) + this.get(2, 2) * m1.get(2, 1);
        lm22 = this.get(2, 0) * m1.get(0, 2) + this.get(2, 1) * m1.get(1, 2) + this.get(2, 2) * m1.get(2, 2);

        this.set(0, 0, lm00);
        this.set(0, 1, lm01);
        this.set(0, 2, lm02);
        this.set(1, 0, lm10);
        this.set(1, 1, lm11);
        this.set(1, 2, lm12);
        this.set(2, 0, lm20);
        this.set(2, 1, lm21);
        this.set(2, 2, lm22);
    }

    /**
     * Sets this matrix to all zeros.
     */
    public final void setZero() {
        set(0, 0, 0.0);
        set(0, 1, 0.0);
        set(0, 2, 0.0);

        set(1, 0, 0.0);
        set(1, 1, 0.0);
        set(1, 2, 0.0);

        set(2, 0, 0.0);
        set(2, 1, 0.0);
        set(2, 2, 0.0);

    }

    /**
     * Negates the value of this matrix: this = -this.
     */
    public final void negate() {
        this.set(0, 0, -this.get(0, 0));
        this.set(0, 1, -this.get(0, 1));
        this.set(0, 2, -this.get(0, 2));

        this.set(1, 0, -this.get(1, 0));
        this.set(1, 1, -this.get(1, 1));
        this.set(1, 2, -this.get(1, 2));

        this.set(2, 0, -this.get(2, 0));
        this.set(2, 1, -this.get(2, 1));
        this.set(2, 2, -this.get(2, 2));
    }

    /**
     * Invert this matrix
     * @return this if successful, null otherwise
     */
    public ArrayBasedMatrix3d invert() {
        return invert(this);
    }

    /**
     * Invert the source matrix and put the result into the destination matrix
     * @param destOrg The destination matrix, or null if a new one is to be created
     * @return The inverted matrix if successful, null otherwise
     */
    public ArrayBasedMatrix3d invert(ArrayBasedMatrix3d destOrg) {
        double determinant = this.determinant();

        if (determinant != 0) {
            ArrayBasedMatrix3d dest = destOrg;
            if (dest == null) {
                dest = new ArrayBasedMatrix3d();
            }

            /* do it the ordinary way
             *
             * inv(A) = 1/det(A) * adj(T), where adj(T) = transpose(Conjugate Matrix)
             *
             * m00 m01 m02
             * m10 m11 m12
             * m20 m21 m22
             */
            double determinantInv = 1f / determinant;

            // get the conjugate matrix
            double t00 = this.get(1, 1) * this.get(2, 2) - this.get(1, 2) * this.get(2, 1);
            double t01 = -this.get(1, 0) * this.get(2, 2) + this.get(1, 2) * this.get(2, 0);
            double t02 = this.get(1, 0) * this.get(2, 1) - this.get(1, 1) * this.get(2, 0);
            double t10 = -this.get(0, 1) * this.get(2, 2) + this.get(0, 2) * this.get(2, 1);
            double t11 = this.get(0, 0) * this.get(2, 2) - this.get(0, 2) * this.get(2, 0);
            double t12 = -this.get(0, 0) * this.get(2, 1) + this.get(0, 1) * this.get(2, 0);
            double t20 = this.get(0, 1) * this.get(1, 2) - this.get(0, 2) * this.get(1, 1);
            double t21 = -this.get(0, 0) * this.get(1, 2) + this.get(0, 2) * this.get(1, 0);
            double t22 = this.get(0, 0) * this.get(1, 1) - this.get(0, 1) * this.get(1, 0);

            dest.set(0, 0, t00 * determinantInv);
            dest.set(1, 1, t11 * determinantInv);
            dest.set(2, 2, t22 * determinantInv);
            dest.set(0, 1, t10 * determinantInv);
            dest.set(1, 0, t01 * determinantInv);
            dest.set(2, 0, t02 * determinantInv);
            dest.set(0, 2, t20 * determinantInv);
            dest.set(1, 2, t21 * determinantInv);
            dest.set(2, 1, t12 * determinantInv);
            return dest;
        } else {
            throw new IllegalStateException("matrix is not invertible");
        }
    }

    /**
     * Multiply this matrix by the tuple t and place the result
     * back into the tuple (t = this*t).
     * @param t  the tuple to be multiplied by this matrix and then replaced
     */
    public final void transformed(Vector3d t) {
        double x;
        double y;
        double z;
        x = get(0, 0) * t.getX() + get(0, 1) * t.getY() + get(0, 2) * t.getZ();
        y = get(1, 0) * t.getX() + get(1, 1) * t.getY() + get(1, 2) * t.getZ();
        z = get(2, 0) * t.getX() + get(2, 1) * t.getY() + get(2, 2) * t.getZ();
        t.set(x, y, z);
    }

    /**
     * Multiply this matrix by the tuple t and and place the result
     * into the tuple "result" (result = this*t).
     * @param t  the tuple to be multiplied by this matrix
     * @return the tuple into which the product is placed
     */
    public final Vector3d createTransformed(BaseVector3d t) {
        double x;
        double y;
        x = getM00() * t.getX() + getM01() * t.getY() + getM02() * t.getZ();
        y = getM10() * t.getX() + getM11() * t.getY() + getM12() * t.getZ();
        return new Vector3d(x, y, getM20() * t.getX() + getM21() * t.getY() + getM22() * t.getZ());
    }

    @Override
    public final double getM00() {
        return get(0, 0);
    }

    /**
     * Set the first matrix element in the first row.
     *
     * @param m00 The m00 to set.
     */
    public final void setM00(double m00) {
        this.set(0, 0, m00);
    }

    @Override
    public final double getM01() {
        return get(0, 1);
    }

    /**
     * Set the second matrix element in the first row.
     *
     * @param m01 The m01 to set.
     */
    public final void setM01(double m01) {
        this.set(0, 1, m01);
    }

    @Override
    public final double getM02() {
        return get(0, 2);
    }

    /**
     * Set the third matrix element in the first row.
     *
     * @param m02 The m02 to set.
     */
    public final void setM02(double m02) {
        this.set(0, 2, m02);
    }

    @Override
    public final double getM10() {
        return get(1, 0);
    }

    /**
     * Set first matrix element in the second row.
     *
     * @param m10 The m10 to set.
     */
    public final void setM10(double m10) {
        this.set(1, 0, m10);
    }

    @Override
    public final double getM11() {
        return get(1, 1);
    }

    /**
     * Set the second matrix element in the second row.
     *
     * @param m11 The m11 to set.
     */
    public final void setM11(double m11) {
        this.set(1, 1, m11);
    }

    @Override
    public final double getM12() {
        return get(1, 2);
    }

    /**
     * Set the third matrix element in the second row.
     *
     * @param m12 The m12 to set.
     */
    public final void setM12(double m12) {
        this.set(1, 2, m12);
    }

    @Override
    public final double getM20() {
        return get(2, 0);
    }

    /**
     * Set the first matrix element in the third row.
     *
     * @param m20 The m20 to set.
     */
    public final void setM20(double m20) {
        this.set(2, 0, m20);
    }

    @Override
    public final double getM21() {
        return get(2, 1);
    }

    /**
     * Set the second matrix element in the third row.
     *
     * @param m21 The m21 to set.
     */
    public final void setM21(double m21) {
        this.set(2, 1, m21);
    }

    @Override
    public final double getM22() {
        return get(2, 2);
    }

    /**
     * Set the third matrix element in the third row.
     *
     * @param m22 The m22 to set.
     */
    public final void setM22(double m22) {
        this.set(2, 2, m22);
    }

}
