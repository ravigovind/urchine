###############################################################################
#                      cannibal genetic algorithm (in r)                      #
###############################################################################
# Copyright (c) 2008-2007 Jonathan Lucas Reddinger <lucas@wingedleopard.net>  #
#                                                                             #
# Permission to use, copy, modify, and distribute this software for any       #
# purpose with or without fee is hereby granted, provided that the above      #
# copyright notice and this permission notice appear in all copies.           #
#                                                                             #
# THE SOFTWARE IS PROVIDED "AS IS" AND THE AUTHOR DISCLAIMS ALL WARRANTIES    #
# WITH REGARD TO THIS SOFTWARE INCLUDING ALL IMPLIED WARRANTIES OF            #
# MERCHANTABILITY AND FITNESS. IN NO EVENT SHALL THE AUTHOR BE LIABLE FOR     #
# ANY SPECIAL, DIRECT, INDIRECT, OR CONSEQUENTIAL DAMAGES OR ANY DAMAGES      #
# WHATSOEVER RESULTING FROM LOSS OF USE, DATA OR PROFITS, WHETHER IN AN       #
# ACTION OF CONTRACT, NEGLIGENCE OR OTHER TORTIOUS ACTION, ARISING OUT OF     #
# OR IN CONNECTION WITH THE USE OR PERFORMANCE OF THIS SOFTWARE.              #
###############################################################################

rm( list = ls(all = TRUE) )         # clear all objects

set.seed(37984732)

# set parameters

num_gens = 100                      ##### !!! look what happens if this is 500
pop_size = 30                       # must be an even number!
num_bits = 8
prob_mut = .01

###############################################################################
# define the functions
###############################################################################

initialize = function() {

    # this function returns a matrix with dimensions [num_bits, pop_size]
    # to look at individual 3, look at the column where j=3
    # this resultant vector can be evaluated as base 2 digits

    pop_matrix = matrix(  rbinom( (pop_size * num_bits) , 1, .5),
                          num_bits,
                          pop_size
                       )

    pop_matrix
}



evaluate = function(gene) {

    # this function accepts vectors of bits (called genes)
    # and evaluates them, returning the fitness

    (bit_values = 2^(rev(0:(num_bits-1))))
    value = sum(gene * bit_values)

    #####################
    # FITNESS FN DOMAIN #
    #####################
    min = 0
    max = 2*pi
   
    # this value will range from 0 to 2^num_bits
    # so lets scale and shift this to the desired range
    x = (value / 2^num_bits) * (max - min) + min

    ####################
    # FITNESS FUNCTION #
    ####################
    eval = -cos(x)

    eval = round(eval, 6)
}



select = function() {

    # window pop_fitness values
    win_fitness = pop_fitness - min(pop_fitness)

    # normalize pop_fitness so that we get a round-roulette
    # probability for each individual to be selected
    prob_fitness = win_fitness / sum(win_fitness)

    # rank() can make win_fitness linear, then can transform with sqrt() &c.

    # sort the population by decreasing prob_fitness
    pop_matrix = pop_matrix[,order(-prob_fitness)]

    # sort the corresponding prob_fitness vector
    prob_fitness = sort(prob_fitness, decreasing = TRUE)

    pop_winners = vector(mode = "integer", length = pop_size)

    # pop_winners = sample(1:pop_size, pop_size, prob_fitness, replace=TRUE)
    # 

    for (winner_slot in 1:pop_size) {

        cutoff = runif(1, 0, 1)
        cum_prob = 0
        candidate = 1

        for (candidate in 1:pop_size) {

            cum_prob = cum_prob + prob_fitness[candidate]

            if (cutoff <= cum_prob) {

                pop_winners[winner_slot] = candidate
                break

            }

        }

        # use sample(n, 1:n, prob, replace=TRUE)

    }

    # i don't want to sort the winners. if i do,
    # then crossover isn't very interesting at all.
    # pop_winners = sort(pop_winners)

    winner_matrix = matrix(0, num_bits, pop_size)

    for (winner in 1:pop_size) {

        #for (bit in 1:num_bits) {

            #winner_matrix[bit, winner] = pop_matrix[bit, pop_winners[winner]]
            winner_matrix[, winner] = pop_matrix[, pop_winners[winner]]

        #}

    }

    winner_matrix

}



crossover = function() {

    child_matrix = matrix(0, num_bits, pop_size)

    for (couple in 1:(pop_size/2)) {

        cut_point = sample(1:num_bits, 1)

        for (bit in 1:num_bits) {

            if (bit <= cut_point) {
                child_matrix[bit, 2*couple-1] = pop_matrix[bit, 2*couple-1]
                child_matrix[bit, 2*couple-0] = pop_matrix[bit, 2*couple-0]
            } else {
                child_matrix[bit, 2*couple-1] = pop_matrix[bit, 2*couple-0]
                child_matrix[bit, 2*couple-0] = pop_matrix[bit, 2*couple-1]

            # use vector[a:b] notation

            }

        }

    }

    child_matrix

}



mutate = function() {

    #       only need this...
    #    mutmatrix = matrix(rbinom( ... use mut_rt ...)...)
    #    newmatrix = ((pop_matrix + mutmatrix) %% 2)

    for (column in 1:pop_size) {

        for (bit in 1:num_bits) {

            if (runif(1, 0, 1) < prob_mut) {

                if (pop_matrix[bit, column] == 1) {

                    pop_matrix[bit, column] = 0

                } else {

                    pop_matrix[bit, column] = 1

                }

            }

        }

    }

    pop_matrix

}

###############################################################################
# run the functions on the population
###############################################################################

pop_max = vector(mode = "numeric", length = num_gens)
pop_avg = vector(mode = "numeric", length = num_gens)

# generate a population matrix [num_bits, pop_size]
pop_matrix = initialize()

for (gen in 1:num_gens) {

    # evaluate the genes (columns) as vectors
    # returns a fitness vector [1, pop_size]
    pop_fitness = apply(pop_matrix, 2, 'evaluate')

    # or can use this --
    # pop_values = t(bit_values) %*% pop_matrix

    # log some statistics
    pop_max[gen] = max(pop_fitness)
    pop_avg[gen] = round(mean(pop_fitness), 6)

    # replace pop_matrix with the parents
    pop_matrix = select()

    # replace pop_matrix with the children
    pop_matrix = crossover()

    # replace pop_matrix with the mutants
    pop_matrix = mutate()

}

plot(1:num_gens, c(rep(0,num_gens-2),min(pop_avg),max(pop_max)), type="n")
points(1:num_gens, pop_max[1:num_gens], col="red")
points(1:num_gens, pop_avg[1:num_gens], col="blue")

