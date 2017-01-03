/*******************************************************************************
 *
 * TRIQS: a Toolbox for Research in Interacting Quantum Systems
 *
 * Copyright (C) 2016 I. Krivenko
 *
 * TRIQS is free software: you can redistribute it and/or modify it under the
 * terms of the GNU General Public License as published by the Free Software
 * Foundation, either version 3 of the License, or (at your option) any later
 * version.
 *
 * TRIQS is distributed in the hope that it will be useful, but WITHOUT ANY
 * WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
 * FOR A PARTICULAR PURPOSE. See the GNU General Public License for more
 * details.
 *
 * You should have received a copy of the GNU General Public License along with
 * TRIQS. If not, see <http://www.gnu.org/licenses/>.
 *
 ******************************************************************************/
#include <triqs/utility/first_include.hpp>
#include <signal.h>
#include <triqs/utility/exceptions.hpp>

#include "signal_handler.hpp"

namespace realevol {
namespace signal_handler {

bool initialized = false;

void callback(int signal) {
 stop();
 if(signal == SIGINT) TRIQS_KEYBOARD_INTERRUPT;
 else                 TRIQS_RUNTIME_ERROR;
}

void start() {
 if (initialized) return;
 static struct sigaction action;
 memset(&action, 0, sizeof(action));
 action.sa_handler = callback;
 sigaction(SIGINT, &action, NULL);
 sigaction(SIGTERM, &action, NULL);
 sigaction(SIGXCPU, &action, NULL);
 sigaction(SIGQUIT, &action, NULL);
 sigaction(SIGUSR1, &action, NULL);
 sigaction(SIGUSR2, &action, NULL);
 sigaction(SIGSTOP, &action, NULL);
 initialized = true;
}

void stop() {
 initialized = false;
}

}}
