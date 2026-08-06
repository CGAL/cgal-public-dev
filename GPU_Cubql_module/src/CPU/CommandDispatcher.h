#pragma once

#include <string>
#include <sstream>
#include <unordered_map>
#include <functional>
#include <iostream>
#include <vector>
#include <iomanip>

class CommandDispatcher {
public:
    using CommandHandler = std::function<void(std::istringstream&)>;

    // Registers a new command with explicit usage syntax and description
    void registerCommand(const std::string& name, 
                         const std::string& usageSyntax, 
                         const std::string& helpDesc, 
                         CommandHandler handler) {
        commands[name] = {usageSyntax, helpDesc, handler};
        orderedNames.push_back(name);
    }

    // Maps an alias to an existing command
    void registerAlias(const std::string& alias, const std::string& targetCommand) {
        aliases[alias] = targetCommand;
    }

    // Parses raw CLI input line and routes to registered handler
    void execute(const std::string& inputLine) {
        std::istringstream iss(inputLine);
        std::string cmd;
        
        if (!(iss >> cmd)) {
            return; // Ignore empty lines
        }

        // Check for aliases first
        auto aliasIt = aliases.find(cmd);
        if (aliasIt != aliases.end()) {
            cmd = aliasIt->second;
        }

        auto it = commands.find(cmd);
        if (it != commands.end()) {
            it->second.handler(iss);
        } else {
            std::cout << "Unknown command: '" << cmd << "'. Type 'help' for a list of commands.\n";
        }
    }

    // Prints formatted command help with usage parameters
    void printHelp() const {
        std::cout << "\n======================================================================\n";
        std::cout << "                         AVAILABLE COMMANDS                           \n";
        std::cout << "======================================================================\n";
        for (const auto& name : orderedNames) {
            const auto& data = commands.at(name);
            std::cout << "  " << std::left << std::setw(12) << name;
            if (!data.usageSyntax.empty()) {
                std::cout << " " << data.usageSyntax;
            }
            std::cout << "\n      - " << data.helpDesc << "\n\n";
        }
        std::cout << "======================================================================\n\n";
    }

private:
    struct CommandData {
        std::string usageSyntax;
        std::string helpDesc;
        CommandHandler handler;
    };

    std::unordered_map<std::string, CommandData> commands;
    std::unordered_map<std::string, std::string> aliases;
    std::vector<std::string> orderedNames;
};