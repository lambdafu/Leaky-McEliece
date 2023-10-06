#ifndef STOCKFISH_LEAKFILE_H
#define STOCKFISH_LEAKFILE_H

#include <iostream>
#include <stdexcept>
#include <string>
#include <fstream>

static int hexdigit(char c) {
    if (c >= '0' && c <= '9') {
        return c - '0';
    } else if (c >= 'a' && c <= 'f') {
        return c - 'a' + 10;
    } else if (c >= 'A' && c <= 'F') {
        return c - 'A' + 10;
    }
    return 0;
}

struct Leakfile {
    std::string filename;

    // The KEM mechanism from which the leak comes.
    std::string kem;

    // The leak bits in ASCII.
    std::string leak;

    // The public key.
    std::string pk;

    // The secret key (for verification).
    std::string sk;

    void parseLine(std::string line, std::string& key, std::string& value) {
        if (line.empty())
            return;
        if (line[0] == '#')
            return;

        // Remove carriage return character if present
        line.erase(line.find_last_not_of(" \t\r") + 1);

        std::size_t equalPos = line.find('=');
        if (equalPos != std::string::npos) {
            key = line.substr(0, equalPos);
            value = line.substr(equalPos + 1);
            // Remove leading/trailing whitespace from key and value
            key.erase(0, key.find_first_not_of(" \t"));
            key.erase(key.find_last_not_of(" \t") + 1);
            value.erase(0, value.find_first_not_of(" \t"));
            value.erase(value.find_last_not_of(" \t") + 1);
            // Remove quotes from value if present
            if (value.size() >= 2 && value.front() == '"' && value.back() == '"') {
                value = value.substr(1, value.size() - 2);
            } else if (value.size() >= 2 && value.front() == '\'' && value.back() == '\'') {
                value = value.substr(1, value.size() - 2);
            }
        }
    }


    Leakfile(const std::string& fname) {
        filename = fname;

        std::ifstream file(fname);
        if (!file)
            throw std::runtime_error("Error opening file: " + fname);

        std::string line;
        while (std::getline(file, line)) {
            std::string key, value;
            parseLine(line, key, value);
            // Keys: kem, sk, pk, ct, ss, leak
            if (key == "kem") {
                kem = value;
            } else if (key == "leak") {
                leak = value;
            } else if (key == "pk") {
                pk = value;
            } else if (key == "sk") {
                sk = value;
            }
        }
        file.close();
    }
};
#endif //STOCKFISH_LEAKFILE_H
