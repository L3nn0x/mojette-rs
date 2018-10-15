use super::std;

pub trait Support<T> {
    fn get_data(&self, l: usize, k: usize) -> &T;
    fn get_data_mut(&mut self, l: usize, k: usize) -> &mut T;
    fn set_data(&mut self, l: usize, k: usize, t: T);
    fn height(&self) -> usize;
    fn width(&self) -> usize;
    fn set_height(&mut self, height: usize);
    fn set_width(&mut self, width: usize);
}

pub struct RectangularSupport<T: Default + std::clone::Clone> {
    height: usize,
    width: usize,
    pub data: Vec<T>
}

pub struct LinearSupport<T: Default + std::clone::Clone> {
    pub data: Vec<T>
}

impl<T> RectangularSupport<T> where T: Default + std::clone::Clone {
    pub fn new(height: usize, width: usize) -> RectangularSupport<T> {
        RectangularSupport {
            height: height,
            width: width,
            data: vec![T::default(); height * width]
        }
    }
}

impl<T> Support<T> for RectangularSupport<T> where T: Default + std::clone::Clone {
    fn get_data(&self, l: usize, k: usize) -> &T {
        &self.data[l * self.width + k]
    }

    fn get_data_mut(&mut self, l: usize, k: usize) -> &mut T {
        &mut self.data[l * self.width + k]
    }

    fn set_data(&mut self, l: usize, k: usize, t: T) {
        self.data[l * self.width + k] = t;
    }

    fn height(&self) -> usize {
        self.height
    }

    fn width(&self) -> usize {
        self.width
    }

    fn set_height(&mut self, height: usize) {
        self.height = height;
        self.data.resize(self.height * self.width, T::default());
    }

    fn set_width(&mut self, width: usize) {
        self.width = width;
        self.data.resize(self.height * self.width, T::default());
    }
}

impl<T> LinearSupport<T> where T: Default + std::clone::Clone {
    pub fn new(size: usize) -> LinearSupport<T> {
        LinearSupport {
            data: vec![T::default(); size]
        }
    }
}

impl<T> Support<T> for LinearSupport<T> where T: Default + std::clone::Clone {
    fn get_data(&self, _: usize, k: usize) -> &T {
        &self.data[k]
    }

    fn get_data_mut(&mut self, _: usize, k: usize) -> &mut T {
        &mut self.data[k]
    }

    fn set_data(&mut self, _: usize, k: usize, t: T) {
        self.data[k] = t;
    }

    fn height(&self) -> usize {
        1
    }

    fn width(&self) -> usize {
        self.data.len()
    }

    fn set_height(&mut self, _: usize) {
    }

    fn set_width(&mut self, width: usize) {
        self.data.resize(width, T::default());
    }
}