#ifndef FRACTION_H
#define FRACTION_H


class error : public std::exception {
public:
    explicit error() = default;
};

class fraction_error : public error {
public:
    explicit fraction_error() = default;
};

class matrix_error : public error {
public:
    explicit matrix_error() = default;
};

class resistive_network_error : public error {
public:
    explicit resistive_network_error() = default;
};

class fraction {
private:
    long long numerator, denominator;
    bool sign;

    long long gcd(long long lhs, long long rhs) {
        if (lhs > rhs) { std::swap(lhs, rhs); }
        if (lhs == 0) {
            return rhs;
        } else {
            return gcd(rhs % lhs, lhs);
        }
    }

    void reduction() {
        if (numerator == 0) {
            denominator = 1;
            return;
        }
        long long tmp = gcd(numerator, denominator);
        numerator /= tmp;
        denominator /= tmp;
    }

public:

    fraction() {
        numerator = 0, denominator = 1;
        sign = true;
    }

    fraction(long long num) {
        numerator = std::abs(num);
        denominator = 1;
        sign = (num >= 0);
    }

    fraction(long long numerator_, long long denominator_) {
        if (denominator_ == 0) {
            throw fraction_error();
        }
        numerator = std::abs(numerator_);
        denominator = std::abs(denominator_);
        sign = !((numerator_ >= 0) ^ (denominator_ >= 0));
        reduction();
    }

    fraction(long long numerator_, long long denominator_, bool sign_) {
        if (denominator_ == 0) {
            throw fraction_error();
        }
        numerator = numerator_;
        denominator = denominator_;
        sign = sign_;
        reduction();
    }

    friend std::istream &operator>>(std::istream &in, fraction &obj) {
        std::string data;
        long long numerator_, denominator_, pos = -1;
        in >> data;
        for (int i = 0; i < data.length(); ++i) {
            if (data[i] == '/') {
                pos = i;
                break;
            }
        }
        try {
            if (pos != -1) {
                numerator_ = std::stoi(data.substr(0, pos));
                denominator_ = std::stoi(data.substr(pos + 1, data.length() - pos - 1));
            } else {
                numerator_ = std::stoi(data);
                denominator_ = 1;
            }
        } catch (...) {
            throw fraction_error();
        }
        if (denominator_ == 0) {
            throw fraction_error();
        }
        obj.numerator = std::abs(numerator_);
        obj.denominator = std::abs(denominator_);
        obj.sign = !((numerator_ >= 0) ^ (denominator_ >= 0));
        obj.reduction();
        return in;
    }

    friend std::ostream &operator<<(std::ostream &out, const fraction &obj) {
        if (!obj.sign && obj.numerator != 0) {
            out << '-';
        }
        if (obj.denominator == 1 || obj.numerator == 0) {
            out << obj.numerator;
        } else {
            out << obj.numerator << '/' << obj.denominator;
        }
        return out;
    }

    friend fraction operator+(const fraction &lhs, const fraction &rhs) {
        if (!(lhs.sign ^ rhs.sign)) {
            return {lhs.numerator * rhs.denominator + lhs.denominator * rhs.numerator,
                    lhs.denominator * rhs.denominator, lhs.sign};
        } else if (lhs.sign) {
            return {lhs.numerator * rhs.denominator - lhs.denominator * rhs.numerator,
                    lhs.denominator * rhs.denominator};
        } else {
            return {-lhs.numerator * rhs.denominator + lhs.denominator * rhs.numerator,
                    lhs.denominator * rhs.denominator};
        }
    }

    friend fraction operator-(const fraction &lhs, const fraction &rhs) {
        if (lhs.sign ^ rhs.sign) {
            if (lhs.sign) {
                return {lhs.numerator * rhs.denominator + lhs.denominator * rhs.numerator,
                        lhs.denominator * rhs.denominator, true};
            } else {
                return {lhs.numerator * rhs.denominator + lhs.denominator * rhs.numerator,
                        lhs.denominator * rhs.denominator, false};
            }
        } else if (lhs.sign) {
            return {lhs.numerator * rhs.denominator - lhs.denominator * rhs.numerator,
                    lhs.denominator * rhs.denominator};
        } else {
            return {-lhs.numerator * rhs.denominator + lhs.denominator * rhs.numerator,
                    lhs.denominator * rhs.denominator};
        }
    }

    friend fraction operator*(const fraction &lhs, const fraction &rhs) {
        return {lhs.numerator * rhs.numerator, lhs.denominator * rhs.denominator,
                !(lhs.sign ^ rhs.sign)};
    }

    friend fraction operator/(const fraction &lhs, const fraction &rhs) {
        if (rhs.numerator == 0) {
            throw fraction_error();
        }
        return {lhs.numerator * rhs.denominator, lhs.denominator * rhs.numerator,
                !(lhs.sign ^ rhs.sign)};
    }

    fraction &operator=(const fraction &obj) = default;

    bool operator==(const fraction &obj) const {
        if (numerator == 0 && obj.numerator == 0) {
            return true;
        }
        return !(sign ^ obj.sign) && numerator == obj.numerator && denominator == obj.denominator;
    }
};

#endif //FRACTION_H

