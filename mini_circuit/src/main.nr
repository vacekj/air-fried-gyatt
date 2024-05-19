fn main(x: Field, y: pub Field) {
    assert(x * 4 == y * 3);
    assert(x != y);

    skibidy(x, y);
}

fn skibidy(x: Field, y: Field) {
    assert(y - x == 2);
}
