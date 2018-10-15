use super::std;

fn gcd<T>(t1: T, t2: T) -> T 
    where T:
        std::cmp::PartialEq +
        std::ops::Rem +
        std::convert::From<<T as std::ops::Rem>::Output> +
        std::convert::From<i8> +
        Copy
{
    if t2 != T::from(0) {
        gcd(t2, T::from(t1 % t2))
    } else {
        t1
    }
}

pub struct Projection<T: Default + std::clone::Clone> {
    p: i16,
    q: u16,
    pub bins: Vec<T>
}

pub type Transform<T> = Vec<Projection<T>>;

impl<T> Projection<T> where T: Default + std::clone::Clone {
    pub fn new(p: i16, q: u16, size: usize) -> Projection<T> {
        let mut proj = Projection {
            p: p,
            q: q,
            bins: vec![T::default(); size]
        };
        proj.reduce();
        proj
    }

    fn reduce(&mut self) {
        let g = gcd(self.p, self.q as i16);
        self.p /= g;
        self.q /= g as u16;
    }

    pub fn p(&self) -> i16 {
        self.p
    }

    pub fn q(&self) -> u16 {
        self.q
    }

    pub fn set_p(&mut self, p: i16) {
        self.p = p;
        self.reduce();
    }

    pub fn set_q(&mut self, q: u16) {
        self.q = q;
        self.reduce();
    }
}

impl<T> Default for Projection<T> where T: Default + std::clone::Clone {
    fn default() -> Projection<T> {
        Projection {
            p: 0,
            q: 1,
            bins: Vec::new()
        }
    }
}

impl<T> std::clone::Clone for Projection<T> where T: Default + std::clone::Clone {
    fn clone(&self) -> Projection<T> {
        Projection {
            p: self.p,
            q: self.q,
            bins: self.bins.clone()
        }
    }
}