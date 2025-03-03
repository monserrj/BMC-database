def make_tables():
    print("I am making tables")


def create_db():
    print("I am creating a database")
    make_tables()


def read_filetype1():
    print("I am reading filetype1")


def read_filetype2():
    print("I am reading filetype2")


def load_data():
    print("I am loading data")
    read_filetype1()
    read_filetype2()


def main():
    create_db()
    load_data()


if __name__ == "__main__":
    main()
