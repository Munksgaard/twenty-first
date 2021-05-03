#![allow(clippy::many_single_char_names)]
use std::fmt;
use std::ops::Add;
use std::ops::Div;
use std::ops::Mul;
use std::ops::Rem;
use std::ops::Sub;

#[derive(Debug, Clone, PartialEq)]
pub struct PrimeField {
    q: i128,
}

impl PrimeField {
    pub fn new(q: i128) -> Self {
        Self { q }
    }
}

#[derive(Debug, Clone, PartialEq, Copy)]
pub struct PrimeFieldElement<'a> {
    pub value: i128,
    pub field: &'a PrimeField,
}

impl fmt::Display for PrimeFieldElement<'_> {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        write!(f, "{} mod {}", self.value, self.field.q)
    }
}

impl<'a> PrimeFieldElement<'a> {
    pub fn new(value: i128, field: &'a PrimeField) -> Self {
        Self {
            value: (value % field.q + field.q) % field.q,
            field,
        }
    }

    pub fn legendre_symbol(&self) -> i128 {
        self.mod_pow((self.field.q - 1) / 2).value
    }

    fn same_field_check(&self, other: &PrimeFieldElement, operation: &str) {
        if self.field.q != other.field.q {
            panic!(
                "Operation {} is only defined for elements in the same field. Got: q={}, p={}",
                operation, self.field.q, self.field.q
            );
        }
    }

    fn egcd(mut x: i128, mut y: i128) -> (i128, i128, i128) {
        let (mut a0, mut a1, mut b0, mut b1) = (1, 0, 0, 1);

        while y != 0 {
            let (q, r) = (x / y, x % y);
            let (c, d) = (a0 - q * a1, b0 - q * b1);

            x = y;
            y = r;
            a0 = a1;
            a1 = c;
            b0 = b1;
            b1 = d;
        }
        (x, a0, b0)
    }
    //     fn egcd_helper(x: i128, y: i128, mut a: i128, mut b: i128) -> (i128, i128, i128) {
    //         // We do this by subtracting the smaller number from the biggest
    //         // as many times as is needed before we hit zero
    //         let x_is_bigger = x > y;
    //         let (mut max, min) = (
    //             if x_is_bigger { x } else { y },
    //             if x_is_bigger { y } else { x },
    //         );
    //         let mut prev = max;
    //         while max > 0 {
    //             prev = max;
    //             max -= min;
    //         }
    //         if max == 0 {
    //             (prev, a, b)
    //         } else {
    //             egcd_helper(prev, min, a, b)
    //         }
    //     }

    //     egcd_helper(x, y, 0, 0)
    // }

    pub fn inv(&self) -> Self {
        let (_, _, a) = Self::egcd(self.field.q, self.value);
        Self {
            field: self.field,
            value: (a % self.field.q + self.field.q) % self.field.q,
        }
    }

    pub fn mod_pow(&self, pow: i128) -> Self {
        let mut acc = Self {
            value: 1,
            field: self.field,
        };
        let res = self.clone();

        for i in 0..128 {
            acc = acc.clone() * acc.clone();
            let set: bool = pow & (1 << (128 - 1 - i)) != 0;
            if set {
                acc = acc * res.clone();
            }
        }
        acc
    }
}

impl<'a> Add for PrimeFieldElement<'a> {
    type Output = Self;

    fn add(self, other: Self) -> Self {
        self.same_field_check(&other, "add");
        Self {
            value: (self.value + other.value) % self.field.q,
            field: self.field,
        }
    }
}

impl<'a> Sub for PrimeFieldElement<'a> {
    type Output = Self;

    fn sub(self, other: Self) -> Self {
        self.same_field_check(&other, "sub");
        Self {
            value: ((self.value - other.value) % self.field.q + self.field.q) % self.field.q,
            field: self.field,
        }
    }
}

impl<'a> Mul for PrimeFieldElement<'a> {
    type Output = Self;

    fn mul(self, other: Self) -> Self {
        self.same_field_check(&other, "mul");
        Self {
            value: self.value * other.value % self.field.q,
            field: self.field,
        }
    }
}

impl<'a> Div for PrimeFieldElement<'a> {
    type Output = Self;

    fn div(self, other: Self) -> Self {
        self.same_field_check(&other, "div");
        Self {
            value: other.inv().value * self.value % self.field.q,
            field: self.field,
        }
    }
}

impl<'a> Rem for PrimeFieldElement<'a> {
    type Output = Self;

    fn rem(self, other: Self) -> Self {
        self.same_field_check(&other, "rem");
        Self {
            value: 0,
            field: self.field,
        }
    }
}

// p = k*n+1 = 2^32 − 2^20 + 1 = 4293918721
// p-1=2^20*3^2*5*7*13.

#[cfg(test)]
mod test_modular_arithmetic {
    #![allow(clippy::just_underscores_and_digits)]

    #[test]
    fn internal() {
        use super::*;

        println!("{:?}", PrimeFieldElement::egcd(3, 19));

        let field_19 = PrimeField { q: 19 };
        let field_23 = PrimeField { q: 23 };
        let _3_19 = PrimeFieldElement {
            value: 3,
            field: &field_19,
        };
        let _5_19 = PrimeFieldElement {
            value: 5,
            field: &field_19,
        };
        assert_eq!(_3_19.inv(), PrimeFieldElement::new(13, &field_19));
        assert_eq!(
            PrimeFieldElement::new(13, &field_19).mul(PrimeFieldElement {
                value: 5,
                field: &field_19
            }),
            PrimeFieldElement::new(8, &field_19)
        );
        assert_eq!(
            PrimeFieldElement::new(13, &field_19)
                * PrimeFieldElement {
                    value: 5,
                    field: &field_19
                },
            PrimeFieldElement::new(8, &field_19)
        );
        assert_eq!(
            PrimeFieldElement::new(13, &field_19)
                + PrimeFieldElement {
                    value: 5,
                    field: &field_19
                },
            PrimeFieldElement::new(18, &field_19)
        );

        assert_eq!(
            PrimeFieldElement::new(13, &field_19)
                + PrimeFieldElement {
                    value: 13,
                    field: &field_19
                },
            PrimeFieldElement::new(7, &field_19)
        );
        assert_eq!(
            PrimeFieldElement::new(2, &field_23).inv(),
            PrimeFieldElement::new(12, &field_23)
        );
        // let ft_size = 8;
        // let mut impulse_data = vec![ComplexNumber::zero(); ft_size];
        // impulse_data[0] = ComplexNumber::one();
        // let mut impulse = Vector::from(impulse_data);

        // let frequency_domain_dft = dtf_slow(&impulse);

        // // DFT implementation, pulse at one
        // impulse_data = vec![ComplexNumber::zero(); ft_size];
        // impulse_data[0] = ComplexNumber::zero();
        // impulse_data[1] = ComplexNumber::one();
        // let mut impulse_new = Vector::from(impulse_data);
        // let frequency_domain_new_dft = dtf_slow(&impulse_new);

        // impulse_data = vec![ComplexNumber::zero(); ft_size];
        // impulse_data[0] = ComplexNumber::one();
        // impulse = Vector::from(impulse_data);

        // #[allow(unused_variables)] // Ignore warnings since we only are interested in runtime
        // let frequency_domain_fft = fft(impulse);

        // // FFT implementation, pulse at one
        // impulse_data = vec![ComplexNumber::zero(); ft_size];
        // impulse_data[0] = ComplexNumber::zero();
        // impulse_data[1] = ComplexNumber::one();
        // impulse_new = Vector::from(impulse_data);
        // #[allow(unused_variables)]
        // let frequency_domain_new_fft = fft(impulse_new);
        // println!("ft_size = {}", ft_size);
        // println!("dft_height = {}", frequency_domain_dft.height());
        // println!("fft_height = {}", frequency_domain_fft.height());
        // for i in 0..ft_size {
        //     assert!(
        //         (frequency_domain_dft.get(i) - frequency_domain_fft.get(i)).get_real() < 0.0001
        //     );
        //     assert!(
        //         (frequency_domain_dft.get(i) - frequency_domain_fft.get(i)).get_imaginary()
        //             < 0.0001
        //     );
        //     assert!(
        //         (frequency_domain_new_dft.get(i) - frequency_domain_new_fft.get(i)).get_real()
        //             < 0.0001
        //     );
        //     assert!(
        //         (frequency_domain_new_dft.get(i) - frequency_domain_new_fft.get(i)).get_imaginary()
        //             < 0.0001
        //     );
        // }
    }
}
