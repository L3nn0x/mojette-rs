pub mod projection;
pub mod support;
pub mod mt;

use projection::{Projection, Transform};
use support::{LinearSupport, RectangularSupport, Support};

pub fn katz_criterion<T: Default + std::clone::Clone>(projections: &Transform<T>, support: &Support<T>) -> bool {
    let sum_p = projections.iter().fold(0, |a, p| a + p.p().abs() as usize);
    let sum_q = projections.iter().fold(0, |a, p| a + p.q() as usize);

    if sum_p >= support.width() || sum_q >= support.height() {
        true
    } else {
        false
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn katz_test() {
        let support = RectangularSupport::<u8>::new(3, 4);

        let mut transform = Transform::<u8>::new();
        transform.push(Projection::<u8>::new(1, 0, 0));
        transform.push(Projection::<u8>::new(1, 1, 0));
        transform.push(Projection::<u8>::new(0, 1, 0));
        transform.push(Projection::<u8>::new(-1, 1, 0));

        assert_eq!(true, katz_criterion(&transform, &support));

        let transform = transform[..3].to_vec();

        assert_eq!(false, katz_criterion(&transform, &support));
    }

    #[test]
    fn katz_test2() {
        let support = RectangularSupport::<u8>::new(3, 3);

        let mut transform = Transform::<u8>::new();
        transform.push(Projection::<u8>::new(1, 0, 0));
        transform.push(Projection::<u8>::new(1, 1, 0));
        transform.push(Projection::<u8>::new(0, 1, 0));
        transform.push(Projection::<u8>::new(-1, 1, 0));

        assert_eq!(true, katz_criterion(&transform, &support));

        let transform = transform[..3].to_vec();

        assert_eq!(false, katz_criterion(&transform, &support));
    }

    #[test]
    fn katz_test3() {
        let support = LinearSupport::<u8>::new(3);

        let mut transform = Transform::<u8>::new();
        transform.push(Projection::<u8>::new(1, 0, 0));
        transform.push(Projection::<u8>::new(1, 1, 0));
        transform.push(Projection::<u8>::new(0, 1, 0));
        transform.push(Projection::<u8>::new(-1, 1, 0));

        assert_eq!(true, katz_criterion(&transform, &support));

        let transform = transform[..1].to_vec();

        assert_eq!(false, katz_criterion(&transform, &support));
    }

    #[test]
    fn katz_test4() {
        let support = RectangularSupport::<u8>::new(4, 3);

        let mut transform = Transform::<u8>::new();
        transform.push(Projection::<u8>::new(1, 0, 0));
        transform.push(Projection::<u8>::new(1, 1, 0));
        transform.push(Projection::<u8>::new(0, 1, 0));
        transform.push(Projection::<u8>::new(-1, 1, 0));

        assert_eq!(true, katz_criterion(&transform, &support));

        let transform = transform[..3].to_vec();

        assert_eq!(false, katz_criterion(&transform, &support));
    }
}