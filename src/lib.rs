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

pub struct RectangularSupport<T: Default + std::clone::Clone> {
    height: usize,
    width: usize,
    data: Vec<T>
}

impl<T> RectangularSupport<T> where T: Default + std::clone::Clone {
    pub fn new(height: usize, width: usize) -> RectangularSupport<T> {
        RectangularSupport {
            height: height,
            width: width,
            data: vec![T::default(); height * width]
        }
    }

    pub fn get_data(&self, l: usize, k: usize) -> &T {
        &self.data[l * self.width + k]
    }

    pub fn get_data_mut(&mut self, l: usize, k: usize) -> &mut T {
        &mut self.data[l * self.width + k]
    }

    pub fn set_data(&mut self, l: usize, k: usize, t: T) {
        self.data[l * self.width + k] = t;
    }

    pub fn height(&self) -> usize {
        self.height
    }

    pub fn width(&self) -> usize {
        self.width
    }

    pub fn set_height(&mut self, height: usize) {
        self.height = height;
        self.data.resize(self.height * self.width, T::default());
    }

    pub fn set_width(&mut self, width: usize) {
        self.width = width;
        self.data.resize(self.height * self.width, T::default());
    }
}

pub fn direct<T: Default + std::clone::Clone + std::ops::AddAssign>(support: RectangularSupport<T>, projections: &mut Transform<T>) {
    let P = support.width(); // k
    let Q = support.height(); // l

    for proj in projections.iter_mut() {
        let p = proj.p();
        let q = proj.q();

        let size = p.abs() as usize * (Q - 1) + q as usize * (P - 1) + 1;

        proj.bins = vec![T::default(); size];

        let offset = if p < 0 {
            (Q - 1) as isize * p as isize
        } else { 0 };

        for l in 0..Q {
            for k in 0..P {
                proj.bins[(k as isize * q as isize + (l as isize * p as isize - offset)) as usize] += support.get_data(l, k).clone();
            }
        }
    }
}

struct Univoc(usize, usize);

pub fn inverse<T: Default + std::clone::Clone + std::ops::AddAssign + std::ops::SubAssign + std::cmp::PartialEq + std::fmt::Display>(support: &mut RectangularSupport<T>, mut projections: Transform<T>) {
    let P = support.width(); // k
    let Q = support.height(); // l
    let mut offsets = Vec::new();
    let mut unitary_projs = Vec::new();
    let mut dietmar_projs = Vec::new();

    // Build unitary and dietmar projections.
    for proj in projections.iter() {
        let p = proj.p();
        let q = proj.q();
        let size = proj.bins.len();

        let offset = if p < 0 {
            (Q - 1) as isize * p as isize
        } else { 0 };
        offsets.push(offset);
        unitary_projs.push(Projection::<usize>::new(p, q, size));
        dietmar_projs.push(Projection::<usize>::new(p, q, size));
    }

    let mut dietmar = 0;
    for l in 0..Q {
        for k in 0..P {
            for (n, proj) in projections.iter().enumerate() {
                let index = (k as isize* proj.q() as isize + (l as isize * proj.p() as isize - offsets[n])) as usize;
                unitary_projs[n].bins[index] += 1;
                dietmar_projs[n].bins[index] += dietmar;
            }
            dietmar += 1;
        }
    }

    // Find initial univocs bins.
    let mut list = std::collections::VecDeque::new();
    for (n, proj) in projections.iter().enumerate() {
        for i in 0..proj.bins.len() {
            if unitary_projs[n].bins[i] == 1 {
                list.push_back(Univoc(n, i));
            }
        }
    }

    // Reconstruct.
    while list.len() != 0 {
        let univoc = list.pop_front().unwrap();
        let dietmar = dietmar_projs[univoc.0].bins[univoc.1];
        let bin = projections[univoc.0].bins[univoc.1].clone();

        if unitary_projs[univoc.0].bins[univoc.1] == 1 {
            // update support
            let l = dietmar / P;
            let k = dietmar - l * P;
            support.set_data(l, k, bin.clone());

            // update projections
            for (n, proj) in projections.iter_mut().enumerate() {
                let index = (k as isize * proj.q() as isize + (l as isize * proj.p() as isize - offsets[n])) as usize;
                proj.bins[index] -= bin.clone();
                dietmar_projs[n].bins[index] -= dietmar;
                unitary_projs[n].bins[index] -= 1;
                if unitary_projs[n].bins[index] == 1 {
                    list.push_back(Univoc(n, index));
                }
            }
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn basic_test() {
        let mut support = RectangularSupport::<usize>::new(3, 3);
        support.data[0] = 3;
        support.data[1] = 8;
        support.data[2] = 2;
        support.data[3] = 4;
        support.data[4] = 5;
        support.data[5] = 42;
        support.data[6] = 90;
        support.data[7] = 34;
        support.data[8] = 1;

        let mut transform = Transform::<usize>::new();
        transform.push(Projection::<usize>::new(1, 0, 0));
        transform.push(Projection::<usize>::new(1, 1, 0));
        transform.push(Projection::<usize>::new(0, 1, 0));
        transform.push(Projection::<usize>::new(-1, 1, 0));

        direct(support, &mut transform);

        let projection = [
            (1, 0, vec![13, 51, 125]),
            (1, 1, vec![3, 12, 97, 76, 1]),
            (0, 1, vec![97, 47, 45]),
            (-1, 1, vec![90, 38, 9, 50, 2])
        ];

        for (p, t) in projection.iter().zip(transform.clone()) {
            assert_eq!(p.0, t.p);
            assert_eq!(p.1, t.q);
            assert_eq!(p.2, t.bins);
        }

        let mut support = RectangularSupport::<usize>::new(3, 3);

        inverse(&mut support, transform);

        let res = vec![3, 8, 2, 4, 5, 42, 90, 34, 1];
        
        assert_eq!(res, support.data);
    }

    #[test]
    fn rectangular_test() {
        let mut support = RectangularSupport::<isize>::new(3, 4);
        support.data[0] = 3;
        support.data[1] = 8;
        support.data[2] = 2;
        support.data[3] = 4;
        support.data[4] = 5;
        support.data[5] = 42;
        support.data[6] = 90;
        support.data[7] = 34;
        support.data[8] = 1;
        support.data[9] = 55;
        support.data[10] = 23;
        support.data[11] = 99;

        let mut transform = Transform::<isize>::new();
        transform.push(Projection::<isize>::new(1, 0, 0));
        transform.push(Projection::<isize>::new(1, 1, 0));
        transform.push(Projection::<isize>::new(0, 1, 0));
        transform.push(Projection::<isize>::new(-1, 1, 0));

        direct(support, &mut transform);

        let mut support = RectangularSupport::<isize>::new(3, 4);

        inverse(&mut support, transform);

        let res = vec![3, 8, 2, 4, 5, 42, 90, 34, 1, 55, 23, 99];
        
        assert_eq!(res, support.data);
    }

    #[test]
    fn rectangular_test2() {
        let mut support = RectangularSupport::<isize>::new(4, 3);
        support.data[0] = 3;
        support.data[1] = 8;
        support.data[2] = 2;
        support.data[3] = 4;
        support.data[4] = 5;
        support.data[5] = 42;
        support.data[6] = 90;
        support.data[7] = 34;
        support.data[8] = 1;
        support.data[9] = 55;
        support.data[10] = 23;
        support.data[11] = 99;

        let mut transform = Transform::<isize>::new();
        transform.push(Projection::<isize>::new(1, 0, 0));
        transform.push(Projection::<isize>::new(1, 1, 0));
        transform.push(Projection::<isize>::new(0, 1, 0));
        transform.push(Projection::<isize>::new(-1, 1, 0));

        direct(support, &mut transform);

        let mut support = RectangularSupport::<isize>::new(4, 3);

        inverse(&mut support, transform);

        let res = vec![3, 8, 2, 4, 5, 42, 90, 34, 1, 55, 23, 99];
        
        assert_eq!(res, support.data);
    }
}