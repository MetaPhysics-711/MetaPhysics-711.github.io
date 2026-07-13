const cloudinary = require('cloudinary').v2

cloudinary.config({
    cloud_name: 'digumuwth',
    secure: true
})

const url = cloudinary.url('屏幕截图_2025-07-15_182731_aopvfz')

console.log(url)
