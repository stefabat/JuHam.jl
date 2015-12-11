# Define parent type of all models
abstract Model

type TBModel <: Model
	bc::AbstractString
    tij::Function
end

